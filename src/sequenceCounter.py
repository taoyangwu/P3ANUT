import anytree
from anytree.exporter import UniqueDotExporter
import json
import pandas as pd
import re
import numpy as np
import logomaker as lm
# from  sequenceCounter import countJsontFile, countCSVFile
import os
from sklearn.cluster import DBSCAN, OPTICS
import csv


#---------------------------------------------------------#
#Function: validEncodings
#Description: Function to return the valid encoding methods, these need to be hardcoded
#Inputs: None
#Outputs: list - a list of valid encoding methods
#---------------------------------------------------------#
def validEncodings():
    return ["ONEHOT", "BLSOUM", "None"]

#---------------------------------------------------------#
#Function: validMethods
#Description: Function to return the valid counting methods, these need to be hardcoded
#Inputs: None
#Outputs: list - a list of valid methods
#---------------------------------------------------------#
def validMethods():
    return ["direct", "DBSCAN", "OPTICS"]

#---------------------------------------------------------#
#Function: countJsonFile
#Description: An function that allows the user to count the sequences in a json file. This
#           function will parse the json file from the paired assembler and count the sequences
#           based on the tags that are provided. The user can also specify the encoding method
#           and the counting method.
#Inputs: jsonFile - string - the path to the json file
#        parseAmino - boolean - a flag to determine if the amino acid sequences should be counted
#        parseDNA - boolean - a flag to determine if the DNA sequences should be counted
#        basedirectory - string - the base directory to output the files
#        **kwargs - dict - a dictionary of optional arguments
#Outputs: None
#---------------------------------------------------------#
def countJsonFile(jsonFile, parseAmino, parseDNA,
                   basedirectory = "", **kwargs):
    
    proteinTag = kwargs.get("proteinDatatag", "proteinSequence")
    dnaTag = kwargs.get("dnaDatatag", "sequences")
    
    with open(jsonFile, 'r') as stream:
        data = json.load(stream)
     
    sequenceArray = [] 
    
    sample = next(iter(data.values()))
    # print(sample.type())
    
    if(parseAmino):
        
        if(proteinTag not in sample):
            print("Error: No protein tag found")
        else:
            for key, value in data.items():
                if(proteinTag in value):
                    sequenceArray.append(value[proteinTag])
                    
            if(len(sequenceArray) == 0):
                print("Error: No  protein found for tag: ", proteinTag)
                return
        
            parse(sequenceArray, baseDirectory=basedirectory, aminoConversion=True, **kwargs)
        
    if(parseDNA):
        if(dnaTag not in sample):
            print("Error: No DNA tag found")
        else:
            for key, value in data.items():
                if(dnaTag in value):
                    sequenceArray.append(value[dnaTag])
                    
            if(len(sequenceArray) == 0):
                print("Error: No  protein found for tag: ", dnaTag)
        
            parse(sequenceArray, baseDirectory=basedirectory, aminoConversion=False, **kwargs)
        
#---------------------------------------------------------#
#Function: countCSVFile
#Description: An function that allows the user to count the sequences in a csv file. This
#           function will parse the csv file and count the sequences based on the tags that are provided.
#           The user can also specify the encoding method and the counting method.
#Inputs: csvFile - string - the path to the csv file
#        aminoConversion - boolean - a flag to determine wether the sequences are amino acids or DNA
#        basedirectory - string - the base directory to output the files
#        **kwargs - dict - a dictionary of optional arguments
#Outputs: None
def countCSVFile(csvFile,basedirectory, aminoConversion = True, **kwargs):
    sequenceArray = []
    
    with open(csvFile, 'r') as file:
        reader = csv.reader(file)
        
        #Skip the header
        reader.__next__()
        for row in reader:
            sequenceArray.append(row[0])
            
    parse(sequenceArray, aminoConversion,basedirectory, **kwargs)

#---------------------------------------------------------#
#Function: parse
#Description: A function that will parse the sequences and count the sequences based on the encoding and counting method
#Inputs: data - list - a list of sequences
#        aminoConversion - boolean - a flag to determine wether the sequences are amino acids or DNA
#        baseDirectory - string - the base directory to output the files
#        **kwargs - dict - a dictionary of optional arguments
#Outputs: None
#---------------------------------------------------------#
def parse(data, aminoConversion = True, baseDirectory = "", **kwargs):
    
    #Loading Arguments from kwargs that are reponsible for  determmionign the functionality of the counting Method
    method = kwargs.get("countingMethod", "direct")
    encoding = kwargs.get("encodingMethod", "ONEHOT")
    
    #This toggle will output of all the sequence that are not matched by the clustering or counting method
    includeUnMatched = kwargs.get("includeUnMatched", True)
    
    #KWARG that can be used to change the start and ending sequence barcodes, used in the regex statements
    sequenceStart = kwargs.get("sequenceStart", "SHSS")
    sequenceEnd = kwargs.get("sequenceEnd", "GGGS")
    
    #The length of the middle sequence in amino acids
    middleMinLength = kwargs.get("middleMinLength", 5)
    middleMaxLength = kwargs.get("middleMaxLength", 10)
    
    #Inorder to increase the effieciecy of the script, a minimum count can be set to filter out the low count sequences after the counting or clustering operation
    minimumCount = kwargs.get("minimumCount", 1)
    purgedCSV = kwargs.get("purgedCSV", True)
    
    #With the sequences of different lengths, the user can set a threshold to determine if the sequences should be outputed
    differentLengthThreshold = kwargs.get("differentLengthThreshold", 0.01)
    differentLengthBase = kwargs.get("differentLengthBase", True)
    
    #A boolean toggle to create a csv file of the unmatched sequences from the regex statements
    unMatchedRegex = kwargs.get("unMatchedRegex", True)
    
    createTree = kwargs.get("createTree", None)
    createLogo = kwargs.get("createLogo", "t")
    
    maxSequenceLength = kwargs.get("maxSequenceLength", 32)
    maxSequenceLength = maxSequenceLength * 3 if not aminoConversion else maxSequenceLength
    
    #Simple lambda to add the base directory to the putput file Names
    addBaseDirectory = lambda x: os.path.join(baseDirectory, x)
    
    #Add Amino or DNA to the output file name
    filePrefix = "Amino_" if aminoConversion else "DNA_"
    
    #Key is the length of the middle sequence, value is the regex statement
    regexStatements = {}
    #Creating a regex detection statement for when the sequence constainds amino acids
    if(aminoConversion):
        
        
        sequenceStart = re.sub(r"((?<!\\)\*)", "\*",sequenceStart.upper())
        sequenceEnd = re.sub(r"((?<!\\)\*)", "\*",sequenceEnd.upper())
        
        regexStatement = '(?:.*)(sequenceStartsequenceMiddlesequenceEnd)'
        regexStatement = regexStatement.replace("sequenceStart", sequenceStart)
        regexStatement = regexStatement.replace("sequenceEnd", sequenceEnd)
        
        for i  in range( middleMaxLength, middleMinLength - 1, -1):
            
            regexStatements[i] = regexStatement.replace("sequenceMiddle", f"(?:[A-Z]{{{i}}})")
    
    else:
        referenceDictionary = {        
            #**dict.fromkeys(['ATA', 'ATC', 'ATT'], 'I'), 
            "I" : "(?:(?:ATA)|(?:ATT)|(?:ATC))",
            # **dict.fromkeys(['ACA', 'ACC', 'ACG', 'ACT'], 'T'), 
            "T" : "(?:(?:ACA)|(?:ACC)|(?:ACG)|(?:ACT))",
            # **dict.fromkeys(['AAC', 'AAT'], 'N'),
            "N" : "(?:(?:AAC)|(?:AAT))",
            # **dict.fromkeys(['AAA', 'AAG'], 'K'), 
            "K" : "(?:(?:AAA)|(?:AAG))",
            # **dict.fromkeys(['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'], 'S'), 
            "S" : "(?:(?:AGC)|(?:AGT)|(?:TCA)|(?:TCC)|(?:TCG)|(?:TCT))",
            # **dict.fromkeys(['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'], 'R'),
            'R' : "(?:(?:AGA)|(?:AGG)|(?:CGA)|(?:CGC)|(?:CGG)|(?:CGT))",
            # **dict.fromkeys(['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'], 'L'), 
            'L' : "(?:(?:CTA)|(?:CTC)|(?:CTG)|(?:CTT)|(?:TTA)|(?:TTG))",            
            # **dict.fromkeys(['CCA', 'CCC', 'CCG', 'CCT'], 'P'),
            'P' : "(?:(?:CCA)|(?:CCC)|(?:CCG)|(?:CCT))",
            # **dict.fromkeys(['CAC', 'CAT'], 'H'), 
            'H' : "(?:(?:CAC)|(?:CAT))",
            # **dict.fromkeys(['CAA', 'CAG', 'TAG'], 'Q'), # In this case, amber codon is set to glutamine
            'Q' : "(?:(?:CAA)|(?:CAG)|(?:TAG))",
            # **dict.fromkeys(['GTA', 'GTC', 'GTG', 'GTT'], 'V'), 
            'V' : "(?:(?:GTA)|(?:GTC)|(?:GTG)|(?:GTT))",
            # **dict.fromkeys(['GCA', 'GCC', 'GCG', 'GCT'], 'A'),
            'A' : "(?:(?:GCA)|(?:GCC)|(?:GCG)|(?:GCT))",
            # **dict.fromkeys(['GAC', 'GAT'], 'D'),
            'D' : "(?:(?:GAC)|(?:GAT))",
            # **dict.fromkeys(['GAA', 'GAG'], 'E'),
            'E' : "(?:(?:GAA)|(?:GAG))",
            # **dict.fromkeys(['GGA', 'GGC', 'GGG', 'GGT'], 'G'),
            'G' : "(?:(?:GGA)|(?:GGC)|(?:GGG)|(?:GGT))",
            # **dict.fromkeys(['TTC', 'TTT'], 'F'),
            'F' : "(?:(?:TTC)|(?:TTT))",
            # **dict.fromkeys(['TAC', 'TAT'], 'Y'),
            'Y' : "(?:(?:TAC)|(?:TAT))",
            # **dict.fromkeys(['TAA', 'TGA'], '*'),
            '*' : "(?:(?:TAA)|(?:TGA))",
            # **dict.fromkeys(['TGC', 'TGT'], 'C'),
            'C' : "(?:(?:TGC)|(?:TGT))",
            # 'TGG':'W', 
            'W' : "(?:(?:TGG))",
            # 'ATG':'M',
            'M' : "(?:(?:ATG))",
        }
        sequenceStart = re.sub(r"((?<!\\)\*)", "\*",sequenceStart.upper())
        sequenceEnd = re.sub(r"((?<!\\)\*)", "\*",sequenceEnd.upper())
        
        regexStatement = '(?:.*)(sequenceStartsequenceMiddlesequenceEnd)'
        
        convertedSequenceStart = ""
        for char in sequenceStart:
            convertedSequenceStart += referenceDictionary.get(char, char)
            
        convertedSequenceEnd = ""
        for char in sequenceEnd:
            convertedSequenceEnd += referenceDictionary.get(char, char)
            
        regexStatement = regexStatement.replace("sequenceStart", convertedSequenceStart)
        regexStatement = regexStatement.replace("sequenceEnd", convertedSequenceEnd)

        for i  in range( middleMaxLength, middleMinLength -1, -1 ):
            modifedStatement = regexStatement.replace("sequenceMiddle", f"(?:[A-Z]{{{i * 3}}})")
            regexStatements[i*3] = modifedStatement
            
    #Regex to split the list into seperate items
    
    if(encoding == "ONEHOT"):
        encoding = "oneHotProtein" if aminoConversion else "oneHotDNA"
        
    #When Using a Direct Method, the sequences are not encoded
    if(method == "direct"):
        encoding = "None-Amino" if aminoConversion else "None-DNA"
    
    #
            
    #Encode the sequences
    print("Encoding Sequences")
    encodedSequences,bases = seqEncoding(encoding, data, maxSequenceLength, **kwargs)
    numOfSequences = len(encodedSequences)
    print("Finished Encoding Sequences")
    
    countDF, unMatched = None, None
    #Count the methods
    if(method == "direct"):
        countDF, unMatched = directComparisionCount(encodedSequences, **kwargs)
    elif(method == "DBSCAN"):
        if encoding == "None" or (encoding == "oneHotProtein" and not aminoConversion) or (encoding == "oneHotDNA" and aminoConversion):
            "Error: Incorrect encoding for DBSCAN"
            encoding = "oneHotProtein" if aminoConversion else "oneHotDNA"
        countDF, unMatched = dbScanCount(encodedSequences, encoding, **kwargs)
    elif(method == "OPTICS"):
        if encoding == "None" or (encoding == "oneHotProtein" and not aminoConversion) or (encoding == "oneHotDNA" and aminoConversion):
            "Error: Incorrect encoding for DBSCAN"
            encoding = "oneHotProtein" if aminoConversion else "oneHotDNA"
        countDF, unMatched = opticsCount(encodedSequences, encoding, **kwargs)
    else:
        print("Error: method not found. Using direct method")
        countDF, unMatched = directComparisionCount(encodedSequences, **kwargs)
        
        
    #Check if the user wants to output the raw comparison data
    if(includeUnMatched and len(unMatched) > 0):
        dataFrameToCSV(unMatched, addBaseDirectory(filePrefix + "unmatched.csv"))

    
    #Purge element that are below the threshold count
    purgedDF = countDF[countDF["m_index"] >= minimumCount]
    
    if(purgedCSV):
        dataFrameToCSV(purgedDF[countDF["m_index"] < minimumCount], addBaseDirectory((filePrefix + "purged.csv")))
        
        
    aboveThreshold = {}
    for key, statement in regexStatements.items():
        
        print(f"Processing {key} length sequences")
        print(f"Regex Statement: {statement}")
        print(f"amount of items in the DF: {len(purgedDF)}")
        #Determine if the statement is in the index (the sequence), output a boolean array
        statementMatch = purgedDF.index.str.contains(statement)
        amountMatched = np.sum(statementMatch)
        print(f"Amount of matched items: {amountMatched}")
        
        if(amountMatched > differentLengthThreshold):
            seq = [re.findall(statement, x)[0] for x in purgedDF[statementMatch].index]
            m_index, s_index = purgedDF[statementMatch].values.T
            
            updatedValues = {}
            for s, m, st in zip(seq, m_index, s_index):
                updatedValues[s] = updatedValues.get(s, 0) + m
                
            newDF = pd.DataFrame.from_dict(updatedValues, orient='index', columns=["m_index"])
            newDF.index.name = "sequence"
            
            newDF.sort_values(by=['m_index'], inplace=True, ascending=False)
            newDF.insert(1, "s_index",0)
                
            
            aboveThreshold[key] = newDF
           
        #Removing the matched sequences from the purgedDF, reduce teh memory and increase the speed of the operations
        purgedDF = purgedDF[~statementMatch]
        print("DF size after purge: ", len(purgedDF))
        
    if(unMatchedRegex):
        dataFrameToCSV(purgedDF, addBaseDirectory(filePrefix + "unmatched.csv"))

    
    for key, value in aboveThreshold.items():
        if(createLogo):
            pass
            createFrequencyDataFrame(value, filePrefix + f"logo_{key}.png", bases, baseDirectory=baseDirectory, **kwargs)
            fileName = filePrefix + f"{key}.csv"
            value.to_csv(fileName)
        if(createTree):
            pass
            createTree(value, **kwargs)
    
    for key, value in aboveThreshold.items():
        print(f"{key} : {len(value)}")
        fileName = filePrefix + f"{key}.csv"
        value.to_csv(addBaseDirectory(fileName))
    
            

    
    

#---------------------------------------------------------#
#Function: directComparisionCount
#Description: A function that will count the sequences based on a direct comparision. This function will count the sequences
#           based on the direct comparision of the sequences. The function will return a dataframe of the sequences and the counts
#Inputs: sequences - list - a list of sequences
#        includeUnMatched - boolean - a flag to determine if the unmatched sequences should be included
#        **kwargs - dict - a dictionary of optional arguments
#Outputs: DataFrame - a dataframe of the sequences and the counts
#---------------------------------------------------------#
def directComparisionCount(sequences, includeUnMatched = False, **kwargs):
    
    sequenceDictionary = {}
    unmatchedDictionary = pd.DataFrame(columns=["sequence", "m_index", "s_index"])
    
    for seq in sequences:
        sequenceDictionary[seq] = sequenceDictionary.get(seq, 0) + 1
    
    sequncesDataFrame = pd.DataFrame.from_dict(sequenceDictionary, orient='index', columns=["m_index"])
    sequncesDataFrame.index.name = "sequence"
    # sequncesDataFrame.reset_index(inplace=True)
    sequncesDataFrame.sort_values(by=['m_index'], inplace=True, ascending=False)
    sequncesDataFrame.insert(1, "s_index",0)
    
    print(unmatchedDictionary)
    

    return sequncesDataFrame, unmatchedDictionary

#---------------------------------------------------------#
#Function: dataFrameToCSV
#Description: A function that will convert a dictionary into a csv file. This is a common function that is used
#           to output the data from the counting methods
#Inputs: inputDictionary - dict - a dictionary of the sequences and the counts
#        fileName - string - the name of the file to output
#        sort - boolean - a flag to determine if the data should be sorted
#Outputs: None
#---------------------------------------------------------#
def dataFrameToCSV(inputDictionary, fileName, sort = False):
    
    if(".csv" in fileName):
            fileName = fileName + ".csv"
    
    #If the input is a dictionary, convert it to a dataframe
    if(isinstance(inputDictionary, pd.DataFrame)):
        unMatchedDataFrame = inputDictionary
    else:
        unMatchedDataFrame = pd.DataFrame.from_dict(inputDictionary, orient='index', columns=["m_index"])
        unMatchedDataFrame.index.name = "sequence"
        unMatchedDataFrame.insert(1, "s_index",0)
        
    
    if(sort):
        unMatchedDataFrame.sort_values(by=['count'], inplace=True, ascending=False)
        
    
        
    unMatchedDataFrame.to_csv(fileName)
    

#---------------------------------------------------------#
#Function: dbScanCount
#Description: A function that will count the sequences based on a DBSCAN clustering method. This function will count the sequences
#           based on the DBSCAN clustering of the sequences. The function will return a dataframe of the sequences and the counts
#Inputs: sequences - list - a list of sequences
#        encoding - string - the encoding method
#        **kwargs - dict - a dictionary of optional arguments
#Outputs: DataFrame - a dataframe of the sequences and the counts
#---------------------------------------------------------#
def dbScanCount(sequences, encoding, **kwargs):
    
    eps = kwargs.get("eps", 0.5)
    min_samples = kwargs.get("min_samples", 5)
        
    print("DBSCAN clustering Starting")
    db = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1).fit(sequences)
    print("DBSCAN clustering Finished")
    
    unmatched = {}
    matched = {}
    coreSamples = {}
    
    print("Starting to decode")
    for i in range(len(db.labels_)):
        label = db.labels_[i]

        if label == -1:
            oringinalSequence = seqDecoding(encoding, sequences[i])
            unmatched[oringinalSequence] = unmatched.get(oringinalSequence, 0) + 1
        else:
            
            if(label not in coreSamples):
                coreSamples[label] = seqDecoding(encoding,sequences[db.core_sample_indices_[label]])
            
            oringinalSequence = coreSamples[label]
            
            matched[oringinalSequence] = matched.get(oringinalSequence, 0) + 1
    
        
    sequenceDf = pd.DataFrame.from_dict(matched, orient='index', columns=["m_index"])
    sequenceDf.index.name = "sequence"
    sequenceDf.sort_values(by=['m_index'], inplace=True, ascending=False)
    sequenceDf.insert(1, "s_index",0)
    
    unmatchedDf = pd.DataFrame.from_dict(unmatched, orient='index', columns=["m_index"])
    unmatchedDf.index.name = "sequence"
    unmatchedDf.sort_values(by=['m_index'], inplace=True, ascending=False)
    unmatchedDf.insert(1, "s_index",0)
    
    return sequenceDf, unmatchedDf

#---------------------------------------------------------#
#Function: opticsCount
#Description: A function that will count the sequences based on a OPTICS clustering method. This function will count the sequences
#           based on the OPTICS clustering of the sequences. The function will return a dataframe of the sequences and the counts
#Inputs: sequences - list - a list of sequences
#        encoding - string - the encoding method
#        **kwargs - dict - a dictionary of optional arguments
#Outputs: DataFrame - a dataframe of the sequences and the counts
#---------------------------------------------------------#
def opticsCount(sequences, encoding, **kwargs):

    min_samples = kwargs.get("min_samples", 5)
        
    print("OPTICS clustering Starting")
    opt = OPTICS(min_samples=min_samples, n_jobs=-1).fit(sequences)
    print("OPTICS clustering Finished")
    
    unmatched = {}
    matched = {}
    
    print("Starting to decode")
    for i in range(len(opt.labels_)):
        label = opt.labels_[i]

        if opt.labels_[i] == -1:
            oringinalSequence = seqDecoding(encoding, sequences[i])
            unmatched[oringinalSequence] = unmatched.get(oringinalSequence, 0) + 1
        else:
            
            oringinalSequence = seqDecoding(encoding,sequences[opt.ordering_[label]])
            
            matched[oringinalSequence] = matched.get(oringinalSequence, 0) + 1
    
        
    sequenceDf = pd.DataFrame.from_dict(matched, orient='index', columns=["m_index"])
    sequenceDf.index.name = "sequence"
    sequenceDf.sort_values(by=['m_index'], inplace=True, ascending=False)
    sequenceDf.insert(1, "s_index",0)
    
    unmatchedDf = pd.DataFrame.from_dict(unmatched, orient='index', columns=["m_index"])
    unmatchedDf.index.name = "sequence"
    unmatchedDf.sort_values(by=['m_index'], inplace=True, ascending=False)
    unmatchedDf.insert(1, "s_index",0)
    
    return sequenceDf, unmatchedDf
    
 
#---------------------------------------------------------#
#Function: createFrequencyDataFrame
#Description: A function that will create a frequency dataframe from a sequence dataframe. This function will take the sequence
#           dataframe and create a frequency dataframe. The function will output a logo and a csv file of the frequency dataframe
#Inputs: sequenceDF - DataFrame - a dataframe of the sequences
#        logoName - string - the name of the logo file
#        validBases - string - the valid bases for the sequences
#        baseDirectory - string - the base directory to output the files
#        **kwargs - dict - a dictionary of optional arguments
#Outputs: None
#---------------------------------------------------------#   
def createFrequencyDataFrame(sequenceDF, logoName, validBases = "*ACDEFGHIKLMNPQRSTVWY", baseDirectory = "", **kwargs):
    
    addBaseDirectory = lambda x: os.path.join(baseDirectory, x)
    
    csvName = kwargs.get("dataFequencyCSV", True)
    
    #dmslogo_colorground is expanded to account for the * base
    logoColor = "colorblind_safe" if len(validBases) < 10 else {
        'A': '#f76ab4',
        'C': '#ff7f00',
        'D': '#e41a1c',
        'E': '#e41a1c',
        'F': '#84380b',
        'G': '#f76ab4',
        'H': '#3c58e5',
        'I': '#12ab0d',
        'K': '#3c58e5',
        'L': '#12ab0d',
        'M': '#12ab0d',
        'N': '#972aa8',
        'P': '#12ab0d',
        'Q': '#972aa8',
        'R': '#3c58e5',
        'S': '#ff7f00',
        'T': '#ff7f00',
        'V': '#12ab0d',
        'W': '#84380b',
        'Y': '#84380b',
        '*' : '#000000'
    }
    
    longestSequence = 0
    sequenceLength = len(sequenceDF.index[0])
    highestBase = 0
    lowestBase = 512
    matrix = []
    
    #Find the highest acsii value in the valid bases
    for base in validBases:
        highestBase = max(highestBase, ord(base))
        lowestBase = min(lowestBase, ord(base)-1)
        
    #Create a reference array that is the size of the valid bases(one hot encoding)
    referenceArray = np.zeros((highestBase + 1 - lowestBase, len(validBases)), dtype=np.uint32)
    for base in validBases:
        referenceArray[ord(base) - lowestBase][validBases.find(base)] = 1
        
    #Create an array that is a maxs size
    frequencyMatrix = np.zeros((sequenceLength, len(validBases)), dtype=np.uint32)
   
    for index, row in sequenceDF.iterrows():   
        # Convert into a numpy array or ints (UTF-8)
        longestSequence = max(longestSequence, len(index))
        sequnceByteArray = bytearray(index, 'utf-8')
        sequence = np.frombuffer(sequnceByteArray, dtype=np.uint8)
        np.subtract(sequence, lowestBase, out=sequence, dtype=np.uint8)
        
        paddedSequence = np.zeros((sequenceLength), dtype=np.uint8)
        #Padd the array to match the size of the matrix
        if(len(sequence) > sequenceLength):
            paddedSequence = sequence[:sequenceLength]
        else:
            paddedSequence[:len(sequence)] = sequence
        
        # Take the reference array and use it to index the sequence array
        tempMatrix = np.take(referenceArray, paddedSequence, axis=0)
        # print(tempMatrix)
        # count = np.left_shift(row["count"],8)
        count = row["m_index"]
        
        # tempMatrix = np.add(tempMatrix, count, dtype=np.uint32)
        tempMatrix = np.multiply(tempMatrix, count).astype(np.uint32)
        #Add the temp matrix to the main matrix
        frequencyMatrix = np.add(frequencyMatrix, tempMatrix).astype(np.uint32)
        # print(m2)

    #Remove the padding from the matrix
    frequencyMatrix = frequencyMatrix[:longestSequence] 
    
    matrixDataFrame = pd.DataFrame(frequencyMatrix, columns=list(validBases))
    
    
    if(csvName):
        fileName = addBaseDirectory(logoName[:-4] + "_table.csv")
        matrixDataFrame.index.name = "Position"
        matrixDataFrame.to_csv(fileName)    
    
    if(logoName):
        matrixDataFrame = pd.DataFrame(frequencyMatrix, columns=list(validBases))
        logo = lm.Logo(matrixDataFrame, color_scheme=logoColor)
        logo.ax.set_ylabel('Frequency')
        logo.ax.set_xlabel('Position')
        logo.ax.set_title('Amino Acid Frequency')
        #logo.show()
        t = addBaseDirectory(logoName)
        logo.ax.figure.savefig(addBaseDirectory(logoName))
 
#---------------------------------------------------------#
#Function: seqEncoding
#Description: A function that will encode the sequences based on the encoding method. This function will take the sequences
#           and encode them based on the encoding method. The function will return a list of the encoded sequences
#Inputs: type - string - the encoding method
#        sequences - list - a list of sequences
#        maxSequenceLength - int - the maximum sequence length
#        **kwargs - dict - a dictionary of optional arguments
#Outputs: list - a list of the encoded sequences
#---------------------------------------------------------#
def seqEncoding(type, sequneces, maxSequenceLength, **kwargs):
    #If statement to generate the correct encoding matrix
    
    
    encoding, bases = None, None
    if(type == "oneHotProtein"):
        encoding, _ = oneHotProteinEncoding()
    elif(type == "oneHotDNA"):
        pass
    elif(type == "BLSOUM"):
        encoding, bases = blousomEncoding()
    elif(type == "None-Amino"):
        return sequneces, "ACDEFGHIKLMNPQRSTVWY*"
    elif(type == "None-DNA"):
        return sequneces, "TGACN"
    else:
        print("Error: encoding type not found")
        return sequneces, "ACDEFGHIKLMNPQRSTVWY*"
        
    
    #Does it make sense to use a numpy array here?
    encodedSequences = []
    
    maxFoundSequenceLength = 0
    
    for seq in sequneces:
        
        maxFoundSequenceLength = max(maxFoundSequenceLength, len(seq))
        
        seqArray = bytearray(seq.upper(), 'utf-8')
        seqArray = np.frombuffer(seqArray, dtype=np.uint8)
        
        #Fit all of the sequences to the max length
        if(len(seqArray) < maxSequenceLength):
            seqArray = np.pad(seqArray, (0, maxSequenceLength - len(seqArray)), 'constant', constant_values=0)
        elif(len(seqArray) > maxSequenceLength):
            seqArray = seqArray[:maxSequenceLength]
        
        encodedSequence = np.take(encoding, seqArray, axis=0)
        encodedSequences.append(encodedSequence.flatten())
    
    print(f" The maxiumum found sequnce length is {maxFoundSequenceLength}") 
    return encodedSequences, bases

    
#---------------------------------------------------------#
#Function: seqDecoding
#Description: A function that will decode the sequences based on the encoding method. This function will take the sequences
#           and decode them based on the encoding method. The function will return a list of the decoded sequences
#Inputs: type - string - the encoding method
#        sequences - list - a list of sequences
#Outputs: list - a list of the decoded sequences
#---------------------------------------------------------#
def seqDecoding(type, sequneces):
    #If statement to generate the correct decoding matrix
    rawEncoding, bases = None, None
    if(type == "oneHotProtein"):
        rawEncoding, bases = oneHotProteinEncoding()
    elif(type == "oneHotDNA"):
        rawEncoding, bases = oneHotCodonEncoding()
    elif(type == "BLSOUM"):
        rawEncoding, bases = blousomEncoding()
    elif(type == "None"):
        return sequneces, None
    else:
        print("Error: encoding type not found")
        #Exception because the user should not be able to get to this point
        return None
    
    encoding = np.zeros((len(bases), len(bases)), dtype=np.int8)
    for i, base in enumerate(bases):
        encoding[i] = rawEncoding[ord(base)]
    
    decodedSequneces = []
    for seq in sequneces:
        #full = np.full((len(seq), len(bases)), np.expand_dims(seq, axis=2))
        #Make it 3D or have a nested for loop to chech the norm fort each base
        full = np.full((len(encoding), len(bases), len(seq)), seq.T)
        # norm = np.linalg.norm(encoding - full - encoding)

        # minIndex = np.argmin(norm)
        # decodedSequneces.append(np.take(bases, minIndex))
        
        t = ""
        for s in seq:
            full = np.full((len(encoding), len(bases)), s)
            norm = np.linalg.norm(full - encoding, axis=1)
            minIndex = np.argmin(norm)
            t += bases[minIndex]
            
        decodedSequneces.append(t)
        
        #full = np.full((len(seq), len(encoding), len(bases)), seq)
        # fullEncoded = np.full((len(seq), len(encoding), len(bases)), encoding)
        # norm = np.linalg.norm(fullEncoded - seq, axis=1)
        
    return decodedSequneces
        
        
#---------------------------------------------------------#
#Function: oneHotProteinEncoding
#Description: A function that will create a one hot encoding matrix for the amino acids. This function will return the encoding
#           matrix and the valid bases
#Inputs: None
#Outputs: tuple - a tuple of the encoding matrix and the valid bases
#---------------------------------------------------------#
def oneHotProteinEncoding():
    validBases = "*ACDEFGHIKLMNPQRSTVWY"
    
    refenceMatrix = np.zeros((ord("Y") + 1, len(validBases)), dtype=np.int8)
    
    for i, base in enumerate(validBases):
        refenceMatrix[ord(base)][i] = 1
        
    return refenceMatrix, validBases

#---------------------------------------------------------#
#Function: oneHotCodonEncoding
#Description: A function that will create a one hot encoding matrix for the codons. This function will return the encoding
#           matrix and the valid bases
#Inputs: None
#Outputs: tuple - a tuple of the encoding matrix and the valid bases
#---------------------------------------------------------#
def oneHotCodonEncoding():
    validBases = "TAGC"
    
    refenceMatrix = np.zeros((ord("T") + 1, len(validBases)), dtype=np.int8)
    
    for i, base in enumerate(validBases):
        refenceMatrix[ord(base)][i] = 1
        
    return refenceMatrix, validBases

#---------------------------------------------------------#
#Function: blousomEncoding
#Description: A function that will create blousum encoding matrix for the amino acids. This function will return the encoding
#           matrix and the valid bases
#Inputs: None
#Outputs: tuple - a tuple of the encoding matrix and the valid bases
#---------------------------------------------------------#
def blousomEncoding():
    validBases = "CSTAGPDEQNHKRMILVWYF"
    
    '''
        C | S T A G P | D E Q N | H K R | M I L V | W Y F
    C   9 |-1-1 0-3-3 |-3-4-3-3 |-3-3-3 |-1-1-1-1 |-2-2-2
    S  -1 | 4 1 1 0-1 | 0 0 0 1 |-1-1-1 |-1-2-2-2 |-3-2-2
    T  -1 |-1 5 0-2-1 |-1-1-1 0 |-2-1-1 |-1-1-1 0 |-2-2-2
    A   0 | 1 0 4 0-1 |-2-1-1-2 |-2-1-1 |-1-1-1 0 |-3-2-2
    G  -3 | 0-2 0 6-2 |-1-2-2 0 |-2-2-2 |-3-4-4-3 |-2-3-3
    P  -3 |-1-1-1-2 7 |-1-1-1-2 |-2-2-1 |-2-3-3-2 |-4-3-4
    D  -3 | 0-1-2-1-1 | 6 2 0 1 |-1-2-1 |-3-3-4-3 |-4-3-3
    E  -4 | 0-1-1-2-1 | 2 5 2 0 |-1-1-1 |-2-2-3-2 |-3-2-3
    Q  -3 | 0-1-1-2-1 | 0 2 5 0 | 0 1 1 | 0-3-2-2 |-2-1-3
    N  -3 | 1 0-2 0-2 | 1 0 0 6 | 1 0 0 |-2-3-3-3 |-4-2-3
    H  -3 |-1-2-2-2-2 |-1 0 0 1 | 8 0 1 |-3-3-3-3 |-2 2-1
    R  -3 |-1-1-1-2-2 |-2-1 1 0 | 0 5 2 |-1-2-2-2 |-3-2-3
    K  -3 | 0-1-1-2-1 |-1 1 1 0 |-1 2 5 |-1-3-2-2 |-3-2-3
    M  -1 |-1-1-1-3-2 |-3-2 0-2 |-2-1-1 | 5 1 2 1 |-1-1 0
    I  -1 |-2-1-1-4-3 |-3-3-3-3 |-3-3-3 | 1 4 2 3 |-3-1 0
    L  -1 |-2-1-1-4-3 |-4-3-2-3 |-3-2-2 | 2 2 4 1 |-2-1 0
    V  -1 |-2 0 0-3-2 |-3-2-2-3 |-3-3-2 | 1 3 1 4 |-3-1-1
    W  -2 |-3-2-3-2-4 |-4-3-2-4 |-2-3-3 |-1-3-2-3 |11 2 1
    Y  -2 |-2-2-2-3-3 |-3-2-1-2 | 2-2-2 |-1-1-1-1 | 2 7 3
    F  -2 |-2-2-2-3-4 |-3-3-3-3 |-1-3-3 | 0 0 0 -1| 1 3 6
    '''
    lookUpDictionary = {
    "C" : [9, -1,-1, 0, -3, -3, -3, -4, -3, -3, -3, -3, -3, -1, -1, -1, -1, -2, -2, -2],
    "S" : [-1, 4, 1, 1, 0, -1, 0, 0, 0, 1, -1, -1, 0, -1, -2, -2, -2, -3, -2, -2],
    "T" : [-1, 1, 5, 0, -2, -1, -1, -1, -1, 0, -2, -1, -1, -1, -1, -1, 0, -2, -2, -2],
    "A" : [0,  1, 0, 4, 0, -1, -2, -1, -1, -2, -2, -1, -1, -1, -1, -1, 0, -3, -2, -2],
    "G" : [-3, 0,-2, 0, 6, -2, -1, -2, -2, 0, -2, -2, -2, -3, -4, -4, -3, -2, -3, -3],
    "P" : [-3,-1,-1, -1, -2, 7, -1, -1, -1, -2, -2, -2, -1, -2, -3, -3, -2, -4, -3, -4],
    "D" : [-3, 0, -1, -2, -1, -1, 6, 2, 0, 1, -1, -2, -1, -3, -3, -4, -3, -4, -3, -3],
    "E" : [-4, 0, -1, -1, -2, -1, 2, 5, 2, 0, 0, 0, 1, -2, -3, -3, -2, -3, -2, -3],
    "Q" : [-3, 0, -1, -1, -2, -1, 0, 2, 5, 0, 0, 1, 1, 0, -3, -2, -2, -2, -1, -3],
    "N" : [-3, 1, 0, -2, 0, -2, 1, 0, 0, 6, 1, 0, 0, -2, -3, -3, -3, -4, -2, -3],
    "H" : [-3,-1, -2, -2, -2, -2, -1, 0, 0, 1, 8, 0, -1, -1, -3, -3, -3, -2, 2, -1],
    "K" : [-3,-1, -1, -1, -2, -2, -2, 0, 1, 0, 0, 5, 2, -1, -3, -2, -3, -3, -2, -3],
    "R" : [-3, 0, -1, -1, -2, -1, -1, 1, 1, 0, -1, 2, 5, -1, -3, -2, -2, -3, -2, -3],
    "M" : [-1,-1, -1, -1, -3, -2, -3, -2, 0, -2, -1, -1, -1, 5, 1, 2, 1, -1, -1, 0],
    "I" : [-1,-2, -1, -1, -4, -3, -3, -3, -3, -3, -3, -3, -3, 1, 4, 2, 3, -3, -1, 0],
    "L" : [-1,-2, -1, -1, -4, -3, -4, -3, -2, -3, -3, -2, -2, 2, 2, 4, 1, -2, -1, 0],
    "V" : [-1,-2, 0, 0, -3, -2, -3, -2, -2, -3, -3, -3, -2, 1, 3, 1, 4, -3, -1, -1],
    "W" : [-2,-3, -2, -3, -2, -4, -4, -3, -2, -4, -2, -3, -3, -1, -3, -2, -3, 11, 2, 1],
    "Y" : [-2,-2, -2, -2, -3, -3, -3, -2, -1, -2, 2, -2, -2, -1, -1, -1, -1, 2, 7, 3],
    "F" : [-2,-2, -2, -2, -3, -4, -3, -3, -3, -3, -1, -3, -3, 0, 0, 0, -1, 1, 3, 6],
    }
    refenceMatrix = np.zeros((ord("Y") + 1, len(validBases)), dtype=np.int8)
    
    for i, base in enumerate(validBases):
        refenceMatrix[ord(base)] = lookUpDictionary[base]
        
    return refenceMatrix, validBases
 
#---------------------------------------------------------#
#Function: createTree
#Description: A function that will create a tree from a sequence dataframe. NOTE this function is not actively used
#Inputs: sequenceDF - DataFrame - a dataframe of the sequences
#        **kwargs - dict - a dictionary of optional arguments
#Outputs: Node - a tree node
#---------------------------------------------------------# 
def createTree(sequenceDF, **kwargs):
    root = anytree.Node("root")
    for index, row in sequenceDF.iterrows():
        currentNode = root
        for i in range(len(row["cleanedSequence"])):
            char = row["cleanedSequence"][i]
            
            createChild = True
            for child in currentNode.children:
                if(child.name == char):
                    createChild = False
                    currentNode = child
                    currentNode.count += row["count"]
                    break
            if(createChild):
                currentNode = anytree.Node(char, parent=currentNode, count = row["count"], partialSequence=row["cleanedSequence"][:i+1])
            
    #Look for alterinatives, requires the use of graphviz which is a third party library 
    UniqueDotExporter(root)
    # print(anytree.RenderTree(root))
    
    return root


if(__name__ == "__main__"):
    file = "dev_Tools/P1.merge.csv"
    
    runAmino = True
    runDNA = True
    
    output = "dev_Tools"
    countCSVFile(file, output, sequenceStart = "*", sequenceEnd = "*")
