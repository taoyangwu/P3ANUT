import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats as sp
import re
from os.path import join, splitext, basename


'''
Each file has the entry on every other line
'''

def intermediateProcessing(filePath, includeSurrounding=False, targetBaseLength=36, normalizeCount=True, outputDirectory='', outputAmino= True, createLogoplot=True):
    raw_data = []
    
    with open(filePath, 'r') as file:
        lines = file.readlines()
        raw_data = [lines[i].strip() for i in range(len(lines)) if i % 2 == 1]
        
    sequence_Counts = {}
    for entry in raw_data:
        entry_length = len(entry)
        if entry_length not in sequence_Counts:
            sequence_Counts[entry_length] = {}
        
        if entry not in sequence_Counts[entry_length]:
            sequence_Counts[entry_length][entry] = 1
        else:
            sequence_Counts[entry_length][entry] += 1
        
    for key, value in sequence_Counts.items():
        sequence_Counts[key] = dict(sorted(value.items(), key=lambda item: item[1], reverse=True))    
        
    totalCount = 0
    for key, value in sequence_Counts.items():
        sequenceLengthCount = 0
        for seq, count in value.items():
            sequenceLengthCount += count
            totalCount += count
            
        print(f'Length: {key}, Unique Sequences: {len(value)}, Total Sequences: {sequenceLengthCount}')
    
    if(outputAmino):
        aminoDictionary = {}
        for seq, counts in sequence_Counts[targetBaseLength].items():
            amio_seq = aminoConversion(seq)
            aminoDictionary[amio_seq] = aminoDictionary.get(amio_seq, 0) + counts
    else:
        aminoDictionary = sequence_Counts[targetBaseLength].copy()
        
        
    fileBaseName = splitext(basename(filePath))[0]
    if (createLogoplot and False):
        logoplot(aminoDictionary, fileBase=fileBaseName, outputDirectory=outputDirectory)
        
    if(includeSurrounding):
        #Include the sequences that are one base pair shorter
        for seq, counts in sequence_Counts.get(targetBaseLength - 1, {}).items():
            largest_matched = ""
            debug_Scores = []
            for i in range(len(seq) + 1):
                for base in ['A', 'T', 'C', 'G']:
                    modified_seq = seq[:i] + base + seq[i:]
                    amio_seq = aminoConversion(modified_seq)
                    debug_Scores.append((amio_seq, aminoDictionary.get(amio_seq, -1)))
                    if(aminoDictionary.get(amio_seq, 0) > aminoDictionary.get(largest_matched, 0)):
                        largest_matched = amio_seq
            if(largest_matched != ""):
                aminoDictionary[largest_matched] = aminoDictionary.get(largest_matched, 0) + counts
                totalCount += counts
            pass
                
        for seq, counts in sequence_Counts.get(targetBaseLength + 1, {}).items():
            largest_matched = ""
            for i in range(len(seq)):
                modified_seq = seq[:i] + seq[i+1:]
                amio_seq = aminoConversion(modified_seq)
                if(aminoDictionary.get(amio_seq, 0) > aminoDictionary.get(largest_matched, 0)):
                    largest_matched = amio_seq
            if(largest_matched != ""):
                aminoDictionary[largest_matched] = aminoDictionary.get(largest_matched, 0) + counts
                totalCount += counts
        
    sorted_amino = dict(sorted(aminoDictionary.items(), key=lambda item: item[1], reverse=True))
    
    with open(join(outputDirectory, f"{fileBaseName}.csv"), 'w') as f:
        f.write("sequence,mean,std\n")
        f.write(f"NORMALIZED_ONE_COUNT,{1/totalCount if normalizeCount else 1},0\n")
        for seq, count in sorted_amino.items():
            f.write(f"{seq},{count/totalCount if normalizeCount else count},0\n")
                
                
#------------------------------------------------------------------------------#
# Function Name: aminoConversion()
# Description: This function attemps to convert the sequence into a protein
# Inputs: seq - sequence to convert - np.array
#          **kwargs - a dictionary of optional arguments, generally loaded from a yaml
# Outputs: proteinSequence - the protein sequence - str
#------------------------------------------------------------------------------#
def aminoConversion(seqParameter):
  
    
    conversionDicitionary = {
        "TTT":"F", "TTG":"L", "TTC":"F", "TTA":"L", "TGT":"C", "TGG":"W", "TGC":"C", "TGA":"*", "TCT":"S", "TCG":"S", "TCC":"S", "TCA":"S", "TAT":"Y", "TAG":"*", "TAC":"Y", "TAA":"*",
        "GTT":"V", "GTG":"V", "GTC":"V", "GTA":"V", "GGT":"G", "GGG":"G", "GGC":"G", "GGA":"G", "GCT":"A", "GCG":"A", "GCC":"A", "GCA":"A", "GAT":"D", "GAG":"E", "GAC":"D", "GAA":"E",
        "CTT":"L", "CTG":"L", "CTC":"L", "CTA":"L", "CGT":"R", "CGG":"R", "CGC":"R", "CGA":"R", "CCT":"P", "CCG":"P", "CCC":"P", "CCA":"P", "CAT":"H", "CAG":"Q", "CAC":"H", "CAA":"Q",
        "ATT":"I", "ATG":"M", "ATC":"I", "ATA":"I", "AGT":"S", "AGG":"R", "AGC":"S", "AGA":"R", "ACT":"T", "ACG":"T", "ACC":"T", "ACA":"T", "AAT":"N", "AAG":"K", "AAC":"N", "AAA":"K"
    }
    
    return "".join(conversionDicitionary[t] for t in ["".join([seqParameter[j] for j in range(i, i+3)]) for i in range(0, len(seqParameter) - 2, 3) ])