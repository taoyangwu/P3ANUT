import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats as sp
import re
from os.path import join, splitext, basename
import argparse


'''
Each file has the entry on every other line
'''

FLIP_DICT = str.maketrans('ATCG', 'TAGC')

DNA_BASES = ('A', 'T', 'C', 'G')


#------------------------------------------------------------------------------#
# Function Name: readFastaEntries()
# Description: Reads a legacy, score-less sequencing file and returns the raw
#              sequence strings. This is the only part of this module the FASTA
#              block in the node graph calls: it handles the different possible
#              record layouts (dataEnteryLength lines per record) and the
#              optional reverse complement and motif trimming, but performs no
#              amino conversion, counting or plotting.
# Inputs: filePath - path to the file to read - str
#         dataEnteryLength - lines that make up one record (2 for FASTA,
#                            4 for FASTQ) - int
#         flip_sequence - reverse complement each read - bool
#         motrif_presequence_length - bases trimmed from the front - int
#         motrif_postsequence_length - bases trimmed from the end - int
# Outputs: list of sequence strings
#------------------------------------------------------------------------------#
def readFastaEntries(filePath, dataEnteryLength: int = 2, flip_sequence: bool = False,
                     motrif_presequence_length: int = 0, motrif_postsequence_length: int = 0):

    if dataEnteryLength < 2:
        raise ValueError("dataEnteryLength must be at least 2 - the sequence is read "
                         "from the second line of every record")

    with open(filePath, 'r') as file:
        lines = file.readlines()

    #The sequence always sits on the second line of a record
    raw_data = [lines[i].strip() for i in range(len(lines)) if i % dataEnteryLength == 1]

    trimEnd = -motrif_postsequence_length if motrif_postsequence_length else None

    entries = []
    for entry in raw_data:
        #Skip reads too short to survive the trimming
        if len(entry) <= motrif_presequence_length + motrif_postsequence_length:
            continue

        local_entry = entry[::-1].translate(FLIP_DICT) if flip_sequence else entry
        entries.append(local_entry[motrif_presequence_length:trimEnd])

    return entries


#------------------------------------------------------------------------------#
# Function Name: fastaToPipelineData()
# Description: Wraps readFastaEntries so a legacy file enters the pipeline in the
#              same JSON style representation the Paired Assembler produces,
#              letting every downstream block treat both paths identically.
# Inputs: filePath - path to the file to read - str
#         dnaDatatag - JSON tag the sequence is stored under - str
#         **kwargs - forwarded to readFastaEntries
# Outputs: dict of {recordKey: {dnaDatatag: sequence}}
#------------------------------------------------------------------------------#
def fastaToPipelineData(filePath, dnaDatatag: str = "sequences", **kwargs):
    entries = readFastaEntries(filePath, **kwargs)
    fileBaseName = splitext(basename(filePath))[0]

    return {f"{fileBaseName}_{i}": {dnaDatatag: seq} for i, seq in enumerate(entries)}


#------------------------------------------------------------------------------#
# Function Name: countByLength()
# Description: Buckets sequences by their length and counts each distinct
#              sequence within a bucket, highest count first.
# Inputs: entries - the sequences to bucket - iterable of str
# Outputs: dict of {length: {sequence: count}}
#------------------------------------------------------------------------------#
def countByLength(entries):
    sequence_Counts = {}
    for motif_seq in entries:
        sequence_Counts.setdefault(len(motif_seq), {})
        sequence_Counts[len(motif_seq)][motif_seq] = sequence_Counts[len(motif_seq)].get(motif_seq, 0) + 1

    for key, value in sequence_Counts.items():
        sequence_Counts[key] = dict(sorted(value.items(), key=lambda item: item[1], reverse=True))

    return sequence_Counts


#------------------------------------------------------------------------------#
# Function Name: foldSurroundingSequences()
# Description: Folds "surrounding" sequences - those exactly one base off the
#              target length - into the counted result instead of dropping them.
#              Each one is matched to whichever correctly sized sequence it is
#              closest to, meaning the most common sequence reachable by adding
#              or removing a single base. The Sequence Counter block shares this
#              function so both ingestion paths behave identically.
# Inputs: targetDictionary - counts keyed by the converted target length
#                            sequence, modified in place - dict
#         shorterSequences - {sequence: count} one base below the target - dict
#         longerSequences - {sequence: count} one base above the target - dict
#         convert - maps a DNA sequence onto the key space of targetDictionary
# Outputs: targetDictionary, modified in place
#------------------------------------------------------------------------------#
def foldSurroundingSequences(targetDictionary, shorterSequences=None, longerSequences=None,
                             convert=None):

    convert = aminoConversion if convert is None else convert

    #One base short - try inserting every base at every position
    for seq, counts in (shorterSequences or {}).items():
        largest_matched = ""
        for i in range(len(seq) + 1):
            for base in DNA_BASES:
                candidate = convert(seq[:i] + base + seq[i:])
                if targetDictionary.get(candidate, 0) > targetDictionary.get(largest_matched, 0):
                    largest_matched = candidate

        if largest_matched != "":
            targetDictionary[largest_matched] = targetDictionary.get(largest_matched, 0) + counts

    #One base long - try removing every base in turn
    for seq, counts in (longerSequences or {}).items():
        largest_matched = ""
        for i in range(len(seq)):
            candidate = convert(seq[:i] + seq[i + 1:])
            if targetDictionary.get(candidate, 0) > targetDictionary.get(largest_matched, 0):
                largest_matched = candidate

        if largest_matched != "":
            targetDictionary[largest_matched] = targetDictionary.get(largest_matched, 0) + counts

    return targetDictionary


def intermediateProcessing(filePath, includeSurrounding=False, targetBaseLength=36, normalizeCount=True, outputDirectory='',
                           outputAmino= True, createLogoplot=True, motrif_presequence_length : int = 0, motrif_postsequence_length : int = 0,
                           dataEnteryLength : int = 2, flip_sequence : bool = False):

    raw_data = readFastaEntries(filePath, dataEnteryLength=dataEnteryLength, flip_sequence=flip_sequence,
                                motrif_presequence_length=motrif_presequence_length,
                                motrif_postsequence_length=motrif_postsequence_length)

    sequence_Counts = countByLength(raw_data)

    for key, value in sequence_Counts.items():
        sequenceLengthCount = sum(value.values())
        print(f'Length: {key}, Unique Sequences: {len(value)}, Total Sequences: {sequenceLengthCount}')

    if(outputAmino):
        aminoDictionary = {}
        for seq, counts in sequence_Counts[targetBaseLength].items():
            amio_seq = aminoConversion(seq)
            if amio_seq == "":
                continue
            aminoDictionary[amio_seq] = aminoDictionary.get(amio_seq, 0) + counts
    else:
        aminoDictionary = sequence_Counts[targetBaseLength].copy()
        
        
    fileBaseName = splitext(basename(filePath))[0]
    if (createLogoplot):
        logoplot(aminoDictionary, fileBase=fileBaseName, outputDirectory=outputDirectory)
        
    if(includeSurrounding):
        foldSurroundingSequences(
            aminoDictionary,
            shorterSequences=sequence_Counts.get(targetBaseLength - 1, {}),
            longerSequences=sequence_Counts.get(targetBaseLength + 1, {}),
            convert=aminoConversion if outputAmino else (lambda s: s),
        )

    sorted_amino = dict(sorted(aminoDictionary.items(), key=lambda item: item[1], reverse=True))
    
    file_set_size = len(raw_data)

    with open(join(outputDirectory, f"{fileBaseName}_1.csv"), 'w') as f:
        f.write("sequence,m_index,s_index\n")
        f.write(f"NORMALIZED_ONE_COUNT,{1/file_set_size if normalizeCount else 1},0\n")
        for seq, count in sorted_amino.items():
            f.write(f"{seq},{count/file_set_size if normalizeCount else count},0\n")
                
                
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
    
    if "N" in seqParameter:
        return ""
    
    return "".join(conversionDicitionary[t] for t in ["".join([seqParameter[j] for j in range(i, i+3)]) for i in range(0, len(seqParameter) - 2, 3) ])
    
    
def logoplot(sequenceCounts, fileBase="TempLogoPlot", outputDirectory=".",
             validBases="CS*TAGPDEQNHKRMILVWYF", outputPath=None, dpi=None):
    import logomaker as lm

    if not sequenceCounts:
        raise ValueError("No sequences to build a logo plot from")

    validBasedict = {base: idx for idx, base in enumerate(validBases)}

    matrix = np.zeros((len(validBases), max(len(seq) for seq in sequenceCounts.keys())))

    for seq, count in sequenceCounts.items():
        for position, base in enumerate(seq):
            #A base outside the chosen alphabet carries no information here
            if base in validBasedict:
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

    destination = outputPath or join(outputDirectory, f"{fileBase}_logoplot.png")
    logo.ax.figure.savefig(destination, **({"dpi": dpi} if dpi else {}))
    plt.close(logo.ax.figure)

    return destination
    
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


# --------------------------------------------------------------------------- #
#  CLI entry point                                                             #
# --------------------------------------------------------------------------- #

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="FASTA_fileConversion",
        description="Convert a legacy, score-less sequencing file into the counted "
                    "CSV the rest of P3ANUT consumes, or dump the parsed reads as JSON.",
    )
    p.add_argument('--filePath', type=str, required=True,
                   help='Path to the input file')
    p.add_argument('--outputDirectory', type=str, default='output/',
                   help='Output directory for the processed files and graphs')
    p.add_argument('--dataEnteryLength', type=int, default=2,
                   help='Lines per record: 2 for FASTA, 4 for FASTQ (default: 2)')
    p.add_argument('--flip-sequence', dest='flip_sequence', action='store_true',
                   default=False, help='Reverse complement every read as it is parsed')
    p.add_argument('--motrif_presequence_length', type=int, default=0,
                   help='Length of pre-sequence before motif to trim off')
    p.add_argument('--motrif_postsequence_length', type=int, default=0,
                   help='Length of post-sequence after motif to trim off')
    p.add_argument('--targetBaseLength', type=int, default=36,
                   help='Target base length in DNA base pairs, so 12 amino acids = 36')
    p.add_argument('--includeSurrounding', action='store_true', default=False,
                   help='Fold in sequences one base pair shorter or longer')
    p.add_argument('--no-normalizeCount', dest='normalizeCount', action='store_false',
                   default=True, help='Write raw counts instead of normalized counts')
    p.add_argument('--no-createLogoplot', dest='createLogoplot', action='store_false',
                   default=True, help='Skip the logo plot')
    p.add_argument('--no-outputAmino', dest='outputAmino', action='store_false',
                   default=True, help='Count DNA sequences rather than converting to amino acids')
    p.add_argument('--parse-only', dest='parse_only', metavar='OUTPUT.json', default=None,
                   help='Only run the file parsing step - the part the FASTA block uses - '
                        'and write the parsed reads to this JSON file')
    return p


def main():
    args = _build_parser().parse_args()

    if args.parse_only:
        import json
        data = fastaToPipelineData(
            args.filePath,
            dataEnteryLength=args.dataEnteryLength,
            flip_sequence=args.flip_sequence,
            motrif_presequence_length=args.motrif_presequence_length,
            motrif_postsequence_length=args.motrif_postsequence_length,
        )
        with open(args.parse_only, 'w') as f:
            json.dump(data, f, indent=4)
        print(f"Parsed {len(data)} reads to {args.parse_only}")
        return

    intermediateProcessing(
        args.filePath,
        includeSurrounding=args.includeSurrounding,
        targetBaseLength=args.targetBaseLength,
        normalizeCount=args.normalizeCount,
        outputDirectory=args.outputDirectory,
        outputAmino=args.outputAmino,
        createLogoplot=args.createLogoplot,
        motrif_presequence_length=args.motrif_presequence_length,
        motrif_postsequence_length=args.motrif_postsequence_length,
        dataEnteryLength=args.dataEnteryLength,
        flip_sequence=args.flip_sequence,
    )


if __name__ == "__main__":
    main()
