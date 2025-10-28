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
    if (createLogoplot):
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
    
    
def logoplot(sequenceCounts, fileBase="TempLogoPlot", outputDirectory="."):
    import logomaker as lm
    validBases = "CS*TAGPDEQNHKRMILVWYF"
    
    validBasedict = {base: idx for idx, base in enumerate(validBases)}
    
    matrix = np.zeros((len(validBases), max(len(seq) for seq in sequenceCounts.keys())))
    
    for seq, count in sequenceCounts.items():
        for position, base in enumerate(seq):
            matrix[validBasedict[base], position] += count
            
    matrixDataFrame = pd.DataFrame(matrix.T, columns=list(validBases))
    logo = lm.Logo(matrixDataFrame, color_scheme= {
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
    })
    logo.ax.set_ylabel('Frequency')
    logo.ax.set_xlabel('Position')
    logo.ax.set_title('Amino Acid Frequency')
    logo.ax.figure.savefig(join(outputDirectory, f"{fileBase}_logoplot.png"))
    
def comparisionScatter(file1, file2, point = 25, minCount = 0.003):
    data1 = {}
    remaining_1 = {}
    
    with open(file1, 'r') as f:
        next(f)  # Skip header
        next(f)  # Skip NORMALIZED_ONE_COUNT line
        for i, line in enumerate(f):
            seq, mean, std = line.strip().split(',')
            if(float(mean) > 0.0001 and i < 1000):
                data1[seq] = i + 1
                
            remaining_1[seq] = i + 1
            
    data2 = {}
    remaining_2 = {}
    with open(file2, 'r') as f:
        next(f)  # Skip header
        next(f)  # Skip NORMALIZED_ONE_COUNT line
        for i, line in enumerate(f):
            seq, mean, std = line.strip().split(',')
            if(float(mean) > 0.0001 and i < point):
                data2[seq] = i + 1
            
            remaining_2[seq] = i + 1
            
    common_seqs = set(data1.keys()).union(set(data2.keys()))
    common_seqs = common_seqs.intersection(set(remaining_1.keys())).intersection(set(remaining_2.keys()))
    
    x = [remaining_1.get(seq, point) for seq in data2.keys()]
    y = [remaining_2.get(seq, 1000) for seq in data2.keys()]
    d2_keys = list(data2.keys())
    
    
    sequences_Below = [x[i] - y[i] for i in range(len(x))]
    index = np.where(np.array(sequences_Below) > 0)[0]
    print(index)
    
    x1 = np.arange(1, point)
    y1 = x1 * 1 + 0
    
    cutoffPoint = 0
    for i in range(len(y)):
        t = remaining_2[d2_keys[i]]
        t1 = t < minCount
        if remaining_2[d2_keys[i]] < minCount:
            cutoffPoint = i
            break
    
    plt.scatter(x, y)
    plt.scatter(data1["TKRKHPHRRKYR"], data2["TKRKHPHRRKYR"], color='green', s=100, label='TKRKHPHRRKYR')
    plt.axhline(y=remaining_2[d2_keys[10]], color='blue', linestyle='--', label='Cutoff Line')
    plt.plot(x1, y1, color='red', linestyle='--')
    plt.xlabel('P 1 Sequence Ranking (Lower is Better)')
    plt.ylabel('P  Sequence Ranking (Lower is Better)')
    plt.title('Sequence Index Comparison')
    
    #X axis log
    plt.xscale('log')
    
    plt.show()
    
def comparisionScatter2(file1, file2):
    data1 = {}
    
    with open(file1, 'r') as f:
        next(f)  # Skip header
        next(f)  # Skip NORMALIZED_ONE_COUNT line
        for i, line in enumerate(f):
            seq, mean, std = line.strip().split(',')
            if(float(mean) > 0.0001 and i < 100):
                data1[seq] = i + 1
            
    data2 = {}
    with open(file2, 'r') as f:
        next(f)  # Skip header
        next(f)  # Skip NORMALIZED_ONE_COUNT line
        for i, line in enumerate(f):
            seq, mean, std = line.strip().split(',')
            if(float(mean) > 0.0001 and i < 100):
                data2[seq] = i + 1
            
    common_seqs = set(data1.keys()).intersection(set(data2.keys()))
    
    x = [data1[seq] for seq in common_seqs]
    y = [data2[seq] for seq in common_seqs]
    y2 =[(data2[seq] - data1[seq]) / data1[seq] for seq in common_seqs]
    
    sequences_Below = [x[i] - y[i] for i in range(len(x))]
    index = np.where(np.array(sequences_Below) > 0)[0]
    common_seqs_list = list(common_seqs)
    print(index)
    
    x1 = np.arange(1, 100)
    
    plt.scatter(x, y)
    plt.scatter(data1["TKRKHPHRRKYR"], (data2["TKRKHPHRRKYR"] - data1["TKRKHPHRRKYR"]) / data1["TKRKHPHRRKYR"], color='green', s=100, label='TKRKHPHRRKYR')
    #plt.scatter([data1[common_seqs_list[i]] for i in index], [y2[i] for i in index], color='red', s=50, label='Improved in P4')
   
   
    # plt.plot(x1, x1, color='red', linestyle='--')
    plt.xlabel('P 1 Sequence Ranking (Lower is Better)')
    plt.ylabel('Ratio decrease in ranking')
    plt.title('(F4_index - F1_index) / F1_index Comparison')
    plt.show()
    
def overlappingPatched():
    rect1 = plt.Rectangle((1, 1), 20, 2, color='blue', alpha=0.5)
    rect2 = plt.Rectangle((2, 2), 20, 2, color='red', alpha=0.5)
    
    fig, ax = plt.subplots()
    ax.add_patch(rect1)
    ax.add_patch(rect2)
    # plt.xlim(0, 5)
    plt.ylim(0, 5)
    plt.xscale('log')
    plt.show()
            

                
# intermediateProcessing("dev_Tools/P1.merge.fa", includeSurrounding=True, targetBaseLength=36,
#                        outputDirectory='dev_Tools/', outputAmino=True, normalizeCount=False, createLogoplot=True)
# # print(aminoConversion("GCGCAGCGTTAGCATCCTCATGTGCCTAAGTGTCAG"))

comparisionScatter("dev_Tools/P1.merge.csv", "dev_Tools/P3.merge.csv")
# overlappingPatched()

#Flip sequence

            