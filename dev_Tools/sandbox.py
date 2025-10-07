import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats as sp
import re


'''
Each file has the entry on every other line
'''

def intermediateProcessing(filePath):
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
        
    for key, value in sequence_Counts.items():
        totalCount = 0
        for seq, count in value.items():
            totalCount += count
            
        print(f'Length: {key}, Unique Sequences: {len(value)}, Total Sequences: {totalCount}')
    
    aminoDictionary = {}
    for seq, counts in sequence_Counts[36].items():
        amio_seq = aminoConversion(seq)
        aminoDictionary[amio_seq] = aminoDictionary.get(amio_seq, 0) + counts
        
    sorted_amino = dict(sorted(aminoDictionary.items(), key=lambda item: item[1], reverse=True))
    
    with open('dev_Tools/Amino_36.csv', 'w') as f:
        f.write("Sequence,Count,std\n")
        for seq, count in sorted_amino.items():
            f.write(f"{seq},{count},0\n")
                
                
#------------------------------------------------------------------------------#
# Function Name: aminoConversion()
# Description: This function attemps to convert the sequence into a protein
# Inputs: seq - sequence to convert - np.array
#          **kwargs - a dictionary of optional arguments, generally loaded from a yaml
# Outputs: proteinSequence - the protein sequence - str
#------------------------------------------------------------------------------#
def aminoConversion(seqParameter, desiredStartingProteinBase="RVPFYSHS", desiredEndingProteinBase="RVPFYSHS", aminoBaseRange=3):
  
    
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
    
    conversionDicitionary = {
        "TTT":"F", "TTG":"L", "TTC":"F", "TTA":"L", "TGT":"C", "TGG":"W", "TGC":"C", "TGA":"*", "TCT":"S", "TCG":"S", "TCC":"S", "TCA":"S", "TAT":"Y", "TAG":"*", "TAC":"Y", "TAA":"*",
        "GTT":"V", "GTG":"V", "GTC":"V", "GTA":"V", "GGT":"G", "GGG":"G", "GGC":"G", "GGA":"G", "GCT":"A", "GCG":"A", "GCC":"A", "GCA":"A", "GAT":"D", "GAG":"E", "GAC":"D", "GAA":"E",
        "CTT":"L", "CTG":"L", "CTC":"L", "CTA":"L", "CGT":"R", "CGG":"R", "CGC":"R", "CGA":"R", "CCT":"P", "CCG":"P", "CCC":"P", "CCA":"P", "CAT":"H", "CAG":"Q", "CAC":"H", "CAA":"Q",
        "ATT":"I", "ATG":"M", "ATC":"I", "ATA":"I", "AGT":"S", "AGG":"R", "AGC":"S", "AGA":"R", "ACT":"T", "ACG":"T", "ACC":"T", "ACA":"T", "AAT":"N", "AAG":"K", "AAC":"N", "AAA":"K"
    }
    
    return "".join(conversionDicitionary[t] for t in ["".join([seqParameter[j] for j in range(i, i+3)]) for i in range(0, len(seqParameter) - 2, 3) ])
    
   

        
def conversion(entry):
    
    #Set up conversion arrays that will convert the bases into ints. Additionally it provides additional protection against invalid bases
    baseConversionArray = np.full(ord('T') + 1, 0, dtype=np.uint8)
    baseConversionArray[ord('T')] = 1
    baseConversionArray[ord('G')] = 2
    baseConversionArray[ord('C')] = 3
    baseConversionArray[ord('A')] = 4
    
    #Convert the sequence into a numpy array of ints
    sequnceByteArray = bytearray(entry, 'utf-8')
    sequence = np.frombuffer(sequnceByteArray, dtype=np.uint8)
    
    np.take(baseConversionArray, sequence, out=sequence)
    
    return sequence
                
# intermediateProcessing("dev_Tools/P1.merge.fa")

#Flip sequence
t = np.arr
            