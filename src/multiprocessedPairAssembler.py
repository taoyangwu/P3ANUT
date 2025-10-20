import re
import numpy as np
import json
import time
from Levenshtein import editops
import multiprocessing
from functools import partial
import yaml

#Library for profiling
import cProfile, pstats, datetime, os

#---------------------------------------------------------#
#Function: parse
#Description: Main function for the paired assembler. This function will parse the forward and reverse fastq files
#           and merge the sequences together. The merged sequences will be stored in the data dictionary
#Inputs: forwardFile - string - the path to the forward fastq file
#        ReverseFile - string - the path to the reverse fastq file
#        data - dictionary - the dictionary to store the merged sequences
#        globalPatameters - dictionary - the dictionary to store the global parameters, if it is none it will be set to an empty dictionary
#        addToLogFile - boolean - flag to determine if the metadata should be added to the log file
#        **kwargs - dictionary - additional arguments that can be passed into the function
#Outputs: dictionary - the metadata of the run
#---------------------------------------------------------#
def parse(forwardFile, ReverseFile,data = None, globalPatameters = {}, addToLogFile = True, **kwargs):
    
    print(kwargs)
    
    if(data is None) : data = {}
    
    logFile = globalPatameters.get("logFile", "log.csv")
    errorFile = globalPatameters.get("errorFile", "error.json") 
    
    
    metadata = {}
    metadata["forwardFile"] = forwardFile
    metadata["ReverseFile"] = ReverseFile

    
    startTime = time.time() 
    parseFastqFile(forwardFile, data, **kwargs)
    endTime = time.time()
    metadata["forwardTime"] = endTime - startTime
    d1 = len(data)
    
    #
    if ReverseFile:
        startTime = time.time()
        parseFastqFile(ReverseFile, data, flip=True, **kwargs)
        endTime = time.time()
        metadata["ReverseTime"] = endTime - startTime
        d2 = len(data)
    
        metadata["forwardFileCount"] = d1
        metadata["ReverseFileCount"] = d2
        metadata["uniqueCount"] = d2-d1
    else:
        metadata["ReverseTime"] = 0
        metadata["ReverseFileCount"] = d1
        metadata["uniqueCount"] = 0
        d2 = d1
        
    errorData = []
    startTime = time.time()
    multiprocess(data, errorData=errorData, metadata=metadata, **kwargs)
    endTime = time.time()
    metadata["mergeTime"] = endTime - startTime
    metadata["finalCount"] = len(data)
    metadata["totalTime"] = metadata["forwardTime"] + metadata["ReverseTime"] + metadata["mergeTime"]
    metadata["timeOfRun"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    metadata["droppedCount"] = d2 - metadata["finalCount"]
    
    
    print(f"The program has successfully completed in {metadata['totalTime']} seconds")
    
    if(addToLogFile):
        #Col Names for the log csv
        csvColNames = ["timeOfRun", "forwardFile", "ReverseFile", "forwardTime", "ReverseTime", "forwardFileCount", "ReverseFileCount", "uniqueCount", "mergeTime", "finalCount", "totalTime", "droppedCount"]
        
        #Check if the log file exists, if not create it
        if(not os.path.exists(logFile)):
            with open(logFile, "w") as file:
                file.write(",".join([str(x) for x in csvColNames]) + "\n")
        
        #Adding the latest entry to the log file
        with open(logFile, "a") as file:
            file.write(",".join([str(metadata.get(x, "N/A")) for x in csvColNames]) + "\n")
    
    return metadata

    


#---------------------------------------------------------#
#Function: parseFastqFile
#Description: Function to extract and parse the fastq file. The function will extract the sequence and the quality score
#           and store them in the data dictionary. The function will also check for the presence of the N base and
#           will remove the sequence if it is present
#Inputs: filePath - string - the path to the fastq file
#        data - dictionary - the dictionary to store the sequences
#        flip - boolean - flag to determine if the sequence should be flipped
#        **kwargs - dictionary - additional arguments that can be passed into the function
#Outputs: None
#---------------------------------------------------------#
def parseFastqFile(filePath, data, flip = False, **kwargs ):

    #Load the optional arguments from the kwargs
    cull_maxlength = kwargs.get("cull_maxlength",90)
    cull_minlength = kwargs.get("cull_minlength",20)
    includeN = kwargs.get("includeN",False)
    fileSeperator = kwargs.get("fileSeperator",r"\+")
    multiprocessCode = kwargs.get("multiprocess",False)
    processLimit = None if kwargs.get("processLimit", 0) == 0 else kwargs.get("processLimit", None)
    
    #Set up conversion arrays that will convert the bases into ints. Additionally it provides additional protection against invalid bases
    baseConversionArray = np.full(ord('T') + 1, 0, dtype=np.uint8)
    baseConversionArray[ord('T')] = 1
    baseConversionArray[ord('G')] = 2
    baseConversionArray[ord('C')] = 3
    baseConversionArray[ord('A')] = 4
    
    #Seperate conversion array is needed when the sequence is flipped
    flipBaseConversionArray = np.full(ord('T') + 1, 0, dtype=np.uint8)
    flipBaseConversionArray[ord('T')] = 4
    flipBaseConversionArray[ord('G')] = 3
    flipBaseConversionArray[ord('C')] = 2
    flipBaseConversionArray[ord('A')] = 1
    

    #Load and Read the file
    with open(filePath, "r") as file:
        #Example Entry after regex as a tuple
        #('NB501061:163:HVYLLAFX3:1:11101:1980:1063' - Run ID,  
        # '1:N:0:GTATTATCT+CATATCGTT' - Additional Run Information, 
        # 'TGTAGACTATTCTCACTCTTCTTGTCTGGTTCCTCCGCGTCCGACGTGTGGTGGAGGTTCGGTCGACG', - DNA Sequence 
        # 'AAAAAEE<EE<EEEEEEEEEEAEEEE/EEAEEEEA//AA/EAAEEEEEEEEAEEEEEE/EA<EEEEE6' - DNA Quality Score)
        regexExpression = r"@([A-Z0-9.:\-\s]+)(?:\/*)([A-Z0-9:+\/-]+)([CODONS]{cull_minlength,cull_maxlength})(?:\s+SPLIT\s*)([!-I]{cull_minlength,cull_maxlength})".replace("cull_minlength", str(cull_minlength)).replace("cull_maxlength", str(cull_maxlength))
        regexExpression = regexExpression.replace("CODONS", "ATGCN" if includeN else "ATGC")
        
        cleanedFileSeparator = fileSeperator.replace("\\\\", "\\")
        regexExpression = regexExpression.replace("SPLIT", f"[{cleanedFileSeparator}]?")
        entries = re.findall(regexExpression, file.read())
        
    raw = None
    if(multiprocessCode):
        pool = multiprocessing.Pool(processes=processLimit)
        raw = pool.map(partial(conversion, flip=flip), entries)         
    else:
        raw = [conversion(x, flip=flip) for x in entries]                          
    
    for key, qualityScore, sequenceScore in raw:
        if(key is None):
            continue
        if key not in data:
            data[key] = {
                "primarySequence" : qualityScore,
                "primaryScore" : sequenceScore
            }
        elif "secondarySequence" not in data[key]:
            #Since the additional information is already present, we don't need to add it again
            data[key]["secondarySequence"] = qualityScore
            data[key]["secondaryScore"] = sequenceScore
            

#---------------------------------------------------------#
#Function: conversion
#Description: Function to encodes the sequnce and the quality score into a single array
#Inputs: entry - tuple - the tuple that contains the information about the sequence
#        flip - boolean - flag to determine if the sequence should be flipped
#Outputs: tuple - the tuple that contains the information about the sequence
#---------------------------------------------------------#
def conversion(entry, flip = False):
    
    #------------------------------------------------------------------------------#
    #Inorder to compress the represetantions and speed up the programs, the sequence 
    #and the score are combined into a single array. There are 41 possible scores 
    #and 5 possible bases, so we can combine them into a single array where each 
    #element is a single byte. The score is multiplied by 5 (number of bases) and 
    #added to the base to get the final value. 
    #The following is some sample code to show how the compression is preformed and how it can be decompressed
        # score = 32
        # char = 3 (N : 0, T : 1, G : 2, C : 3, A : 4)
        
        # for i in range(0, 5):
        #     for j in range(0, 41):
        #         mixed = j * 5 + i
        # 
        #         print(i,j, mixed % 5, mixed // 5) 
    #------------------------------------------------------------------------------#
    
    #Set up conversion arrays that will convert the bases into ints. Additionally it provides additional protection against invalid bases
    baseConversionArray = np.full(ord('T') + 1, 0, dtype=np.uint8)
    baseConversionArray[ord('T')] = 1
    baseConversionArray[ord('G')] = 2
    baseConversionArray[ord('C')] = 3
    baseConversionArray[ord('A')] = 4
    
    #Seperate conversion array is needed when the sequence is flipped
    flipBaseConversionArray = np.full(ord('T') + 1, 0, dtype=np.uint8)
    flipBaseConversionArray[ord('T')] = 4
    flipBaseConversionArray[ord('G')] = 3
    flipBaseConversionArray[ord('C')] = 2
    flipBaseConversionArray[ord('A')] = 1
    #Prevent edge conditions where the sequence and the score are not the same length
    if(len(entry[2]) != len(entry[3])):
        #mismatchCount += 1
        return None, None, None
    
    #Convert the sequence into a numpy array of ints
    sequnceByteArray = bytearray(entry[2], 'utf-8')
    sequence = np.frombuffer(sequnceByteArray, dtype=np.uint8)
    
    #Convert the scores into a numpy array of ints
    qualityByteArray = bytearray(entry[3], 'utf-8')
    qualityScore = np.frombuffer(qualityByteArray, dtype=np.uint8)
    
    if(flip):
        sequence = np.flip(sequence)
        qualityScore = np.flip(qualityScore)
        np.take(flipBaseConversionArray, sequence, out=sequence)
    else:
        np.take(baseConversionArray, sequence, out=sequence)
    # sequence = codonConversionFunction(sequence).astype(np.uint8)

    
    #assert(np.array_equal(s2, sequence), "The sequence is not the same")

    
    #The quality score is offset by 33 in the fastq file, so we need to subtract it to get the actual score and allow for compression  
    np.subtract(qualityScore, 33, out=qualityScore)
    
    #Force the scores of the N bases to be 0
    NbaseIndex = np.where(sequence == 0)
    np.put(qualityScore, NbaseIndex, 0)
    
    #Sequence score is used when comparing sequences to determine which sequence is more likely to be correct  
    sequenceScore = np.sum(qualityScore)
    
    #Modify the quality score instead of a new array to make slightly more memeory efficent
    np.multiply(qualityScore, 5, out=qualityScore)
    np.add(qualityScore, sequence, out=qualityScore)
    
    return entry[0], qualityScore, sequenceScore


#---------------------------------------------------------#
#Function: mergeInsertOperation
#Description: Function that handles the merging of two sequences that have an insert operation. The function will
#           determine the most likely location of the error and will merge the sequences together
#Inputs: edits - list - the list of edits needed to convert the primary sequence into the secondary sequence
#        primarySequence - numpy array - the primary sequence
#        primaryScore - numpy array - the primary score
#        secondarySequence - numpy array - the secondary sequence
#        secondaryScore - numpy array - the secondary score
#        **kwargs - dictionary - additional arguments that can be passed into the function
#Outputs: numpy array - the merged sequence
#         numpy array - the merged score
#         numpy array - the codon changes
#---------------------------------------------------------#
def mergeInsertOperation(edits, primarySequence, primaryScore, secondarySequence, secondaryScore, **kwargs):
    
    maxScore = 0
    maxSequence = np.zeros_like(primarySequence, dtype=np.uint8)
    codonChanges = None
    radius = kwargs.get("radius",1)
    scoreOffset = kwargs.get("scoreOffset",1)
    
    for operation, sourcePosition, destenationPosition in edits:
            if(operation == 'insert'):
                tempSequence = np.zeros_like(secondarySequence, dtype=np.uint8)
                
                #There exist and edge case where the insert is at the end of the sequence, so we need to check for that
                if(sourcePosition == len(primarySequence)):
                    tempSequence = primarySequence[:-1]
                    scoreOfInsert = primarySequence[-1] // 5
                else:
                    tempSequence[:sourcePosition] = primarySequence[:sourcePosition]
                    tempSequence[sourcePosition:] = primarySequence[sourcePosition + 1:]
                    scoreOfInsert = primarySequence[sourcePosition] // 5            
                
                #The primary score needs to be updated to account for the removal of the base
                mergeSequence, mergedScore, mergeCodonsChanges = mergeSequencesOperation(tempSequence, primaryScore - scoreOfInsert, 
                                                                     secondarySequence, secondaryScore, 
                                                                     radius, scoreOffset)
                
                
                if(mergedScore > maxScore):
                    maxScore = mergedScore
                    maxSequence = mergeSequence
                    codonChanges = mergeCodonsChanges
                    
    
    return maxSequence, maxScore, codonChanges

#---------------------------------------------------------#
#Function: mergeDeleteOperation
#Description: Function that handles the merging of two sequences that have a delete operation. The function will
#           determine the most likely location of the error and will merge the sequences together
#Inputs: edits - list - the list of edits needed to convert the primary sequence into the secondary sequence
#        primarySequence - numpy array - the primary sequence
#        primaryScore - numpy array - the primary score
#        secondarySequence - numpy array - the secondary sequence
#        secondaryScore - numpy array - the secondary score
#        **kwargs - dictionary - additional arguments that can be passed into the function
#Outputs: numpy array - the merged sequence
#         numpy array - the merged score
#         numpy array - the codon changes
#---------------------------------------------------------#
def mergeDeleteOperation(edits, primarySequence, primaryScore, secondarySequence, secondaryScore, **kwargs):
    radius = kwargs.get("radius",1)
    scoreOffset = kwargs.get("scoreOffset",1)
    
    maxScore = 0
    maxSequence = np.zeros_like(primarySequence, dtype=np.uint8)
    codonChanges = None
    
    for operation, sourcePosition, destenationPosition in edits:
            if(operation == 'delete'):
                tempSequence = np.zeros_like(primarySequence, dtype=np.uint8)
                tempSequence[:destenationPosition] = secondarySequence[:destenationPosition]
                tempSequence[destenationPosition + 1:] = secondarySequence[destenationPosition:]
                
                mergeSequence, mergedScore, mergeCodonsChanges = mergeSequencesOperation(primarySequence, primaryScore, 
                                                                     tempSequence, secondaryScore, 
                                                                     radius, scoreOffset)
                     
                if(mergedScore > maxScore):
                    maxScore = mergedScore
                    maxSequence = mergeSequence
                    codonChanges = mergeCodonsChanges
                    
    
    return maxSequence, maxScore, codonChanges


#---------------------------------------------------------#
#Function: mergeSequencesOperation
#Description: Function that handles the merging of two sequences.
#Inputs: primarySequenceParm - numpy array - the primary sequence
#        primaryScoreParm - numpy array - the primary score
#        secondarySequenceParm - numpy array - the secondary sequence
#        secondaryScoreParm - numpy array - the secondary score
#        radius - int - the radius to use for the surronding scores
#        scoreOffset - int - the offset to use for the surronding scores
#Outputs: numpy array - the merged sequence
#         numpy array - the merged score
#         numpy array - the codon changes
#---------------------------------------------------------#
def mergeSequencesOperation(primarySequenceParm, primaryScoreParm, secondarySequenceParm, secondaryScoreParm, radius, scoreOffset):
    
    #-------------------------------------------------------------------------------------------------------------------#
    #Sample code to better understand the process of generating the surronding score for a sequence
    #In the project code, the values are capped at [-1,1], but it is easier to understand with a larger range
    
    # mydata = np.array([1,2,3,4,5,6,7,8,9,10])
    #Radius = 1
    
    #Get a sliding window sum of the data to check the surronding forward scores
    # t1 = np.convolve(mydata,np.ones(radius+1,dtype=int),'valid')
    #Sliding window sum stops before the end of the array, so we need to add the last few elements using a cumulative sum 
    #Multiple flips([::-1]) are need to oreniate the array correctly
    # t2 = np.cumsum(mydata[len(mydata) - radius:len(mydata)][::-1])[::-1]

    # t = np.append(t1,t2)

    # print(mydata)
    # print(t)

    #Same process for the reverse scores as with the forward scores, but the array needs to be in differenet orentations
    # x1 = np.convolve(mydata[::-1],np.ones(radius+1,dtype=int),'valid')
    # x2 = np.cumsum(mydata[0:radius])

    # x = np.append(x2, x1[::-1])
    # print(x)

    #Combine the forward and reverse scores to get the surronding scores
    # final = np.add(x, t)
    # print(final)
    #-------------------------------------------------------------------------------------------------------------------#
    
    #Conditions
    #1. The peptide is the same and the orginal score is higher - Add the scores to represent an increase in confidence
    #2. The peptide is the same and the scores are the same - Add the scores to represent an increase in confidence
    #3. The peptide is the same and the orginal score is lower - Add the scores to represent an increase in confidence
    #4. The peptide is different and the orginal score is higher - Do nothing
    #5. The peptide is different and the scores are the same - check the surronding scores to find the highest, if tie, use primary
    #6. The peptide is different and the orginal score is lower - Replace the orginal score

    primarySequenceScore = primaryScoreParm
    secondarySequenceScore = secondaryScoreParm
    primarySequence = primarySequenceParm
    secondarySequence = secondarySequenceParm
    #Gaurentee that the primary sequence is the one with the highest score
    if(primaryScoreParm < secondaryScoreParm):
        primarySequence, secondarySequence = secondarySequenceParm, primarySequenceParm
        primaryScore, secondaryScore = secondaryScoreParm, primaryScoreParm
        
    (primarySequenceScore, primarySequenceCodons) = np.divmod(primarySequence, 5)
    (secondarySequenceScore, secondarySequenceCodons) = np.divmod(secondarySequence, 5)
    
    
    diffOfSequences = np.sign(np.subtract(primarySequenceScore, secondarySequenceScore, dtype=np.int8))


    mergedScores = np.convolve(diffOfSequences,np.ones((radius*2)+1,dtype=np.int8),'same')

    #Find the indexies of the mathcing codons
    codonMatches = np.equal(primarySequenceCodons, secondarySequenceCodons)
    codonMissMatches = np.logical_not(codonMatches)
    
    #Generate a new array to track the changes and to track changes
    codonChanges = np.zeros_like(codonMissMatches, dtype=np.uint8)
    np.add(codonChanges, codonMissMatches, out=codonChanges)
    
    
    #Since the score is a base 10 representaion of the probabity of the the codon being correct
    #we can add the score to represent the probability of the codon being correct
    finalSequence = np.zeros_like(primarySequence, dtype=np.uint8)
    np.add(finalSequence, primarySequenceScore, where=codonMatches, out=finalSequence)
    np.add(finalSequence, secondarySequenceScore, where=codonMatches, out=finalSequence)
    
    #Generate the encoding of the codons and the scores for matching codons
    #apply minium to prevent overflow
    np.minimum(finalSequence, 40, out=finalSequence)
    np.multiply(finalSequence, 5, out=finalSequence)
    
    #Adding the codon values to the final sequence
    np.add(finalSequence, primarySequenceCodons, where=codonMatches , out=finalSequence)
    

    #Generate supporting Masks that accomidate for sequnce scores offset
    secondarySequenceMajor = np.add(secondarySequenceScore, scoreOffset)
    #Sence using unsinged ints, we can't have negative numbers, so setting all numbers to at least the error range 
    #will prevent any negative numbers from being generated
    secondarySequneceMinor = np.maximum(secondarySequenceScore, scoreOffset)
    np.subtract(secondarySequneceMinor, scoreOffset, out=secondarySequneceMinor)
    
    #Generate the masks to find which codons need to be inserted
    greaterThan = np.greater(primarySequenceScore, secondarySequenceMajor)
    np.logical_and(greaterThan, codonMissMatches, out=greaterThan)
    
    lessThan  = np.less(primarySequenceScore, secondarySequneceMinor)
    np.logical_and(lessThan, codonMissMatches, out=lessThan)
    #No current NOR function in numpy, so creating NOT of OR
    equal = np.logical_not(np.logical_or(greaterThan, np.logical_or(lessThan, codonMatches)))

    #Apply the mask to the final scores and find which are lower than 0 which means that the secondary sequence is better
    equalSelect = np.select([equal], [mergedScores])
    equalSelectLessThan = np.less(equalSelect, 0)
    
    #Short cut way to calculate the greater than mask since we already have the equal mask and the less than mask
    equalSelectGreaterThan = np.logical_xor(equal, equalSelectLessThan)

    #Using the masks, we can find the codons that need to be inserted
    #merge the equal lists to the final sequence
    np.add(finalSequence, primarySequence, where = greaterThan, out=finalSequence)
    np.add(finalSequence, secondarySequence, where = lessThan, out=finalSequence)
    np.add(finalSequence, primarySequence, where = equalSelectGreaterThan, out=finalSequence)
    np.add(finalSequence, secondarySequence, where = equalSelectLessThan, out=finalSequence)

    
    #the final scores and allows the use of debugging to check the final codons
    (finalSequenceScore, finalSequenceCodons) = np.divmod(finalSequence, 5)
    
    
    return (finalSequence, np.sum(finalSequenceScore), codonChanges)


#---------------------------------------------------------#
#Function: multiprocess
#Description: Function to manage the multiprocessed version of the merge operation
#Inputs: data - dictionary - the dictionary to store the sequences
#        debug - boolean - flag to determine if the debug information should be printed
#        errorData - list - the list to store the error data
#        FileName - string - the name of the file to store the error data
#        metadata - dictionary - the dictionary to store the metadata
#        **kwargs - dictionary - additional arguments that can be passed into the function
#Outputs: None
#---------------------------------------------------------#
def multiprocess(data, debug = False, errorData = [], FileName = "mergeMismatchData.json", metadata = {}, **kwargs):
    radius = kwargs.get("radius",1)
    scoreOffset = kwargs.get("scoreOffset",1)
    multiprocessCode = kwargs.get("multiprocess",False)
    processLimit = None if kwargs.get("processLimit", 0) == 0 else kwargs.get("processLimit", 0)

    errorCount = 0
    idsToDrop = set()
    mappedResults = []
    
    # mappedResults = pool.map(partial(mergeController, debug=debug, **kwargs), data.items())
    # mappedResults = [partial(mergeController, radius=radius, scoreOffset=scoreOffset, debug=debug, **kwargs)(x) for x in data.items()]
    if(multiprocessCode):
        pool = multiprocessing.pool.Pool(processes=processLimit)
        mappedResults = pool.map(partial(mergeController, debug=debug, **kwargs), data.items())
    else:
        for key,value in data.items():
            mappedResults.append(mergeController((key,value), radius=radius, scoreOffset=scoreOffset, debug=debug))
    
    metadata["mappedLength"] = len(mappedResults)
    droppedCount = 0

    for d in mappedResults:
        if(d[0] is None):
            droppedCount += 1

            continue
        key, sequenceString, stringScores, proteinSequence, drop, codonChanges, errordata = d
        if(drop or (sequenceString is None) or (stringScores is None)):
            
            idsToDrop.add(key)
            droppedCount += 1
            continue
        if(errorData is not None and debug):
            errorData.append(errorData)
           
        #Clear the dictionary of any existing data 
        data[key] = {
            "sequences" : sequenceString,
            "scores" : stringScores,
            
            "codonChanges" : "".join([str(x) for x in codonChanges]),
        }
        
        if(proteinSequence is not None):
            data[key]["proteinSequence"] = proteinSequence
    pass        
    for id in idsToDrop:
        data.pop(id)
    pass
    metadata["droppedCount"] = droppedCount


    
    
        
        
            

#---------------------------------------------------------#
#Function: mergeController
#Description: Function to manage the merge operation
#Inputs: keyValue - tuple - the tuple that contains the key and the value
#        debug - boolean - flag to determine if the debug information should be printed
#        **kwargs - dictionary - additional arguments that can be passed into the function
#Outputs: tuple - the tuple that contains the information about the sequence
#---------------------------------------------------------#

def mergeController(keyValue, debug = False, **kwargs):
    radius = kwargs.get("radius",1)
    scoreOffset = kwargs.get("scoreOffset",1)
    lengthMisMatch = 0
    #Return a tuple, (key, maxSequence, maxScore, codonChanges, errorCount(None or {""}))
    #Case to check if there is a secondary sequence to speed up performace
    
    test = lambda x,y, z : print("Error", z) if x is None or y is None else None
    
    key, value = keyValue
    if("secondarySequence" not in value and "primarySequence" in value):
        sequenceString, stringScores, proteinSequence, drop = finalize(key, value["primarySequence"], **kwargs)
        # test(sequenceString, stringScores, 1)
        return (key, sequenceString, stringScores, proteinSequence, None, np.zeros_like(value["primarySequence"]), None) 
    
    if("primarySequence" not in value and "secondarySequence" in value):
        sequenceString, stringScores, proteinSequence, drop = finalize(key, value["secondarySequence"], **kwargs)
        # test(sequenceString, stringScores, 2)
        return (key, sequenceString, stringScores, proteinSequence, None, np.zeros_like(value["secondarySequence"]), None) 
    
    if("primaryScore" not in value and "secondaryScore" not in value):
        print("Error: No primary or secondary score")
        return (None, None, None,None, None, None, None)
    
    #Compare the lengths of the sequence to determine if there is a base deletion or insertion
    if(len(value["primarySequence"]) != len(value["secondarySequence"])):
        
        
        #If the sequences don't match additional steps are needed to determine the most likely location of the error
        maxSequence, maxScore, codonChanges, errorData = mergeMisMatchedLengths(key, value["primarySequence"], value["primaryScore"],
                                                                    value["secondarySequence"], value["secondaryScore"], debug, **kwargs)

        sequenceString, stringScores, proteinSequence, drop = finalize(key, maxSequence, **kwargs)
        return (key, sequenceString, stringScores, proteinSequence, drop, codonChanges, errorData)  
        # #Store the results of the merge
        # value["primarySequence"], value["primaryScore"], value["codonChanges"] = maxSequence, maxScore, codonChanges                                                                    
                        
    else:
        maxSequence, maxScore, codonChanges = mergeSequencesOperation(value["primarySequence"], value["primaryScore"], 
                                value["secondarySequence"], value["secondaryScore"],
                                radius, scoreOffset)
        sequenceString, stringScores, proteinSequence, drop = finalize(key, maxSequence, **kwargs)
        # test(sequenceString, stringScores,3)
        return (key, sequenceString, stringScores, proteinSequence, drop, codonChanges, None)   
    
    
        

#---------------------------------------------------------#
#Function: finalize
#Description: Function to finalize the sequence. It drops the sequence if the quality score is too low and attempts to convert the sequence into a protein
#Inputs: key - string - the id of the sequence
#        sequence - numpy array - the sequence to convert
#        **kwargs - dictionary - additional arguments that can be passed into the function
#Outputs: string - the sequence
#         string - the scores
#         string - the protein sequence
#         boolean - flag to determine if the sequence should be dropped
#---------------------------------------------------------#

def finalize(key, sequence, **kwargs):
    #Load the optional arguments from the kwargs
    minQualityScore = kwargs.get("minQualityScore",18)
    minQualityScoreCount = kwargs.get("minQualityScoreCount",10)
    proteinConversion = kwargs.get("proteinConversion",True)
    aminoBaseRange = kwargs.get("aminoBaseRange", 4)
    
    #Set up an array that will be used to convert the sequence into a string
    baseConversionArray = np.full(5, 0, dtype=np.uint8)
    baseConversionArray[0] = ord('N')
    baseConversionArray[1] = ord('T')
    baseConversionArray[2] = ord('G')
    baseConversionArray[3] = ord('C')
    baseConversionArray[4] = ord('A')
    

    scores, seq = np.divmod(sequence, 5)
        
    lowQualityCount = np.count_nonzero(np.less(scores, minQualityScore))
    
    if(lowQualityCount > minQualityScoreCount):
        #Store the ids to remove and remove them after the loop to prevent errors from changing the size of the dictionary
        return None, None, None, True
    else:
        
        #TODO: Come up with a possible solution to handle a N base 
        #Currenlty if an N base is detected, the sequence is not converted to a protein sequence
        #have a subsutute char for unknown protein
        protein = None
        if(proteinConversion and min(seq) > 0):
            protein = aminoConversion(seq, **kwargs)
            
        stringRep = np.take(baseConversionArray, seq).tobytes()
        np.add(scores, 33, out=scores)
        
        # if(np.max(scores) > 73):
        #     print(seq)
        #     pass

        stringScores = scores.tobytes()
        
        # data[id]["finalSequence"] = stringRep.decode("utf-8")
        # data[id]["finalScore"] = stringScores.decode("utf-8") 
        return stringRep.decode("utf-8"), stringScores.decode("utf-8"), protein, False

#------------------------------------------------------------------------------#
# Function Name: aminoConversion()
# Description: This function attemps to convert the sequence into a protein
# Inputs: seq - sequence to convert - np.array
#          **kwargs - a dictionary of optional arguments, generally loaded from a yaml
# Outputs: proteinSequence - the protein sequence - str
#------------------------------------------------------------------------------#
def aminoConversion(seqParameter, **kwargs):
    
    #Check that there isn't any invalid bases, N case
    if(min(seqParameter) <= 0):
        return None
    
    #Load the optional arguments from the kwargs
    aminoBaseRange = kwargs.get("aminoBaseRange", 3)
    desiredStartingProteinBase = kwargs.get("desiredStartingProteinBase", "RVPFYSHS")
    desiredEndingProteinBase = kwargs.get("desiredEndingProteinBase", "RVPFYSHS")
    proteinSequenceDifference = kwargs.get("proteinSequenceDifference", 2)
    
    #Inorder to convert the sequence into a protein, groups of three codon are converted
    #into a base 4 number. This number then corelates to a specific amino acid. 
    dotMultiplier = np.array([4**2, 4**1, 1**0])
    
    #Since there is already a check to determine that there aren't any N cases( 0 in encoding)
    #we can subtract 1 from the sequence to allow for the sequence to be converted into a base 4 number
    #which correspondingIndex takes up less space
    seq = np.subtract(seqParameter, 1, dtype=np.uint8)
    
    #First row comment denotes the sequence, second row denotes the ascii representation of the amino acid
    correspondingIndex = np.array([
        #TTT, TTG, TTC, TTA, TGT, TGG, TGC, TGA, TCT, TCG, TCC, TCA, TAT, TAG, TAC, TAA
        #"F", "L", "F", "L", "C", "W", "C", "*", "S", "S", "S", "S", "Y", "Q", "Y", "*",
        70,   76,  70,  76,  67,  87,  67,  42,  83,  83,  83,  83,  89,  81,  89,  42,
        
        #GTT, GTG, GTC, GTA, GGT, GGG, GGC, GGA, GCT, GCG, GCC, GCA, GAT, GGG, GAC, GGA
        #"V", "V", "V", "V", "G", "G", "G", "G", "A", "A", "A", "A", "D", "E", "D", "E", 
        86,   86,  86,  86,  71,  71,  71,  71,  65,  65,  65,  65,  68,  69,  68,  69,
        
        #CTT, CTG, CTC, CTA, CGT, CGG, CGC, CGA, CCT, CCG, CCC, CCA, CAT, CAG, CAC, CAA
        #"L", "L", "L", "L", "R", "R", "R", "R", "P", "P", "P", "P", "H", "Q", "H", "Q", 
        76,   76,  76,  76,  82,  82,  82,  82,  80,  80,  80,  80,  72,  81,  72,  81, 
        
        #ATT, ATG, ATC, ATA, AGT, AGG, AGC, AGA, ACT, ACG, ACC, ACA, AAT, AAG, AAC, AAA
        #"I", "M", "I", "I", "S", "R", "S", "R", "T", "T", "T", "T", "N", "K", "N", "K", 
        73,   77,  73,  73,  83,  82,  83,  82,  84,  84,  84,  84,  78,  75, 78,  75
    ])
    
    
    aminoBaseOffset = aminoBaseRange // 3
    # a list to store the priority of the sequences
    #The first slot is considered the highes priority
    priority = [{},{},{}]
    #Loop through the range 
    for i in range(aminoBaseOffset, aminoBaseOffset + aminoBaseRange):
        #Caluclate the sub sequence length, start at the ith index and keep as many as possible
        #How to manage the when the range is 0
        trimmedSeq = seq[i: len(seq) - ((len(seq) - i) % 3)]
        trimmedSeq = np.take(correspondingIndex, np.dot(trimmedSeq.reshape(-1,3), dotMultiplier), mode='clip').astype(np.uint8)
        trimmedSeq = trimmedSeq.tobytes().decode("utf-8")
        
        #Calculating the indexies of the searcged bases
        beginSearch = re.search(f"({desiredStartingProteinBase})(?:[a-zA-Z]+)", trimmedSeq)
        endSearch = re.search(f"(?:[a-zA-Z]*)({desiredEndingProteinBase})(?:[a-zA-Z]+)", trimmedSeq)
        beginIndex = beginSearch.end(1) if beginSearch else 0
        endIndex = endSearch.start(1) if endSearch else len(trimmedSeq)
            
        #Conditions
        #The begining is found but not the end - medium priority -  the length become determining factor
        #The end is found but not the begining - medium priority -  the length become determining factor
        #Neither is found - lowest - priority
        #Both are found - Highest priority
        
        seqDictionary = {
                "sequence": trimmedSeq,
                "beginIndex": beginIndex,
                "length": endIndex - beginIndex,
            }
        
        p = 0
        #Highest
        if(beginSearch and endSearch):
            p = 0
        elif(beginSearch or endSearch):
            p = 1
        else:
            p = 2
            
        if(not priority[p]):
            priority[p] = seqDictionary
        else:
            if(seqDictionary["length"] > priority[p]["length"]):
                priority[p] = seqDictionary
    
    #Returns the highest priority sequence and the assoiated dictionary    
    for p in priority:
        if(p):
            return p["sequence"]
    
    
#---------------------------------------------------------#
#Function: mergeMisMatchedLengths
#Description: Function to handle the merging of two sequences that have a length mismatch
#Inputs: sequenceID - string - the id of the sequence
#        primarySequence - numpy array - the primary sequence
#        primaryScore - numpy array - the primary score
#        secondarySequence - numpy array - the secondary sequence
#        secondaryScore - numpy array - the secondary score
#        **kwargs - dictionary - additional arguments that can be passed into the function
#Outputs: numpy array - the merged sequence
#         numpy array - the merged score
#         numpy array - the codon changes
#         dictionary - the error data
#---------------------------------------------------------#   
def mergeMisMatchedLengths(sequenceID, primarySequence, primaryScore, 
                            secondarySequence, secondaryScore, debug = False, **kwargs):
    
    
    
    #Load the optional arguments from the kwargs
    radius = kwargs.get("radius",1)
    scoreOffset = kwargs.get("scoreOffset",1)
    diffOfLength = len(primarySequence) - len(secondarySequence)
    
    targetLength = kwargs.get("targetLength", 66)
    
    #Determine which sequence is longer
    #Create local copies to allow for local modification and not unknownly effecting
    #the variables
    if(diffOfLength < 0):
        localPrimarySequence, localPrimaryScore = secondarySequence, secondaryScore
        localSecondarySequence, localSecondaryScore = primarySequence, primaryScore
        diffOfLength = abs(diffOfLength)
    else:
        localPrimarySequence, localPrimaryScore = primarySequence, primaryScore
        localSecondarySequence, localSecondaryScore = secondarySequence, secondaryScore
        
    #Print information about the error
    print("Error: The length of the sequences are not the same ")
    print("ID: " + sequenceID)
    print("Primary Sequence: " + str(len(localPrimarySequence)))
    print("Secondary Sequence: " + str(len(localSecondarySequence)))
    
    #Logic to compare the sequences and the most likely location of the error
    #Currenty methodology doesn't handle cases where the error is largers than 2
    if(diffOfLength > 1):
        
        return localPrimarySequence, localPrimaryScore, np.zeros_like(localPrimarySequence, dtype=np.uint8), None
    
    #Variables to store the best sequence and score
    maxScore = 0
    maxSequence = 0
    codonChanges = None
    
    #The levenstein operations needed to convert the primary sequence into the secondary sequence
    edits = editops(localPrimarySequence, localSecondarySequence)
    
    #Cases (primary is the longest, there is difference in length)
    #Assumes 1 distance difference
    #Where there are multiple different inserions/deletions loop through each case and find the best one
    #1. Primary seq matches the target length : replace deletions with a N base w/ score of 0;
    #2. Secondary seq matches the target length : delete insertion possible swap primary and secondary
    #3. Primary seq % 3 == 0 : replace deletions with a N base w/ score of 0;
    #4. Secondary seq % 3 == 0 : delete insertion possible swap primary and secondary
    #5. Else (primary % 3 == 2 and secondary % 3 == 1 or primary % 3 == 1 and secondary % 3 == 2):     
    # Check that a sequence matches rhe taret length
    # remove possible confilcts intoi;l they match target length

    #Case 1 and 3
    if(len(localPrimarySequence) == targetLength or len(localPrimarySequence) % 3 == 0):
        caseNumber = 1 if len(localPrimarySequence) == targetLength else 3
        print(f"Case {caseNumber}")
        #Function is put into a method to reduce code duplication becasue case 5 has to call these methods
        maxSequence, maxScore, codonChanges = mergeDeleteOperation(edits, localPrimarySequence, localPrimaryScore, 
                                                              localSecondarySequence, localSecondaryScore, 
                                                               **kwargs)
                    
    #Case 2 and 4
    elif(len(localSecondarySequence) == targetLength or len(localSecondarySequence) % 3 == 0):
        caseNumber = 2 if len(localSecondarySequence) == targetLength else 4
        print(f"Case {caseNumber}")
        maxSequence, maxScore, codonChanges = mergeInsertOperation(edits, localPrimarySequence, localPrimaryScore, 
                                                              localSecondarySequence, localSecondaryScore, 
                                                              **kwargs)
    #Case 5   
    else:
        print("Case 5")
        caseNumber = 5
        maxSequenceInsert, maxScoreInsert, codonChangesInsert = mergeInsertOperation(edits, localPrimarySequence, localPrimaryScore, 
                                                              localSecondarySequence, localSecondaryScore, 
                                                               **kwargs)
        
        maxSequenceDelete, maxScoreDelete, codonChangesDelete = mergeDeleteOperation(edits, localPrimarySequence, localPrimaryScore, 
                                                              localSecondarySequence, localSecondaryScore,
                                                               **kwargs)
        
        #An additional comparision step is needed to determine which sequence is the best from the insert and delete operations
        if(maxScoreInsert > maxScoreDelete):
            maxSequence, maxScore, codonChanges = maxSequenceInsert, maxScoreInsert, codonChangesInsert
        else:
            maxSequence, maxScore, codonChanges = maxSequenceDelete, maxScoreDelete, codonChangesDelete
                  
    errorData = {
                    "id": sequenceID,
                    "primaryLength": len(primarySequence),
                    "secondaryLength": len(secondarySequence),
                    "primarySequence": [int(x) for x in primarySequence],
                    "secondarySequence": [int(x) for x in secondarySequence],
                    "maxSequence": [int(x) for x in maxSequence],
                    "maxScore": int(maxScore),
                    "caseNumber" : caseNumber
                }
     
    return maxSequence, maxScore, codonChanges, errorData if debug else None

def stripParmDictionary(dict):
    newDict = {}
    
    for key, value in dict.items():
        newDict[key] = value["Value"]
    
    return newDict
    

if(__name__ == "__main__"):
    forFile = "data/example-cut-for.fastq"
    revFile = "data/example-cut-rev.fastq"
    
    # configFile = "configV2.yaml"
    # with open(configFile, 'r') as stream:
    #     try:
    #         yamlSettings = yaml.safe_load(stream)
    #     except yaml.YAMLError as exc:
    #         print(exc)
    #         exit()
    
    # profileParseFastQ("data/PID-1309-GAL-BSA-2-PC_S108_R1_001.fastq")
    data = {}
    meta = parse(forFile, revFile, data=data, multiprocess=False)
    print(meta)
    # print(data)
