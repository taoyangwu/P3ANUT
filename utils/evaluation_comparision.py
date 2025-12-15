import pandas as pd
import csv
import numpy as np
import json
import os

import sys, pathlib
# add project root (two levels up) to sys.path so `from src...` works
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from src.multiprocessedPairAssembler import parseFastqFile, finalize, aminoConversion, parse

def aminoConversion_simplified(seqParameter):
   #GGTGGAGGTTCG
    
    conversionDicitionary = {
        "TTT":"F", "TTG":"L", "TTC":"F", "TTA":"L", "TGT":"C", "TGG":"W", "TGC":"C", "TGA":"*", "TCT":"S", "TCG":"S", "TCC":"S", "TCA":"S", "TAT":"Y", "TAG":"*", "TAC":"Y", "TAA":"*",
        "GTT":"V", "GTG":"V", "GTC":"V", "GTA":"V", "GGT":"G", "GGG":"G", "GGC":"G", "GGA":"G", "GCT":"A", "GCG":"A", "GCC":"A", "GCA":"A", "GAT":"D", "GAG":"E", "GAC":"D", "GAA":"E",
        "CTT":"L", "CTG":"L", "CTC":"L", "CTA":"L", "CGT":"R", "CGG":"R", "CGC":"R", "CGA":"R", "CCT":"P", "CCG":"P", "CCC":"P", "CCA":"P", "CAT":"H", "CAG":"Q", "CAC":"H", "CAA":"Q",
        "ATT":"I", "ATG":"M", "ATC":"I", "ATA":"I", "AGT":"S", "AGG":"R", "AGC":"S", "AGA":"R", "ACT":"T", "ACG":"T", "ACC":"T", "ACA":"T", "AAT":"N", "AAG":"K", "AAC":"N", "AAA":"K"
    }
    
    if "N" in seqParameter:
        return ""
    
    return "".join(conversionDicitionary[t] for t in ["".join([seqParameter[j] for j in range(i, i+3)]) for i in range(0, len(seqParameter) - 2, 3) ])


def evalutate_fastq_file(file, flip = False, start_barcode = "AAA", end_barcode = "TTT", **kwargs):
    
    
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
    
    file_data = {}
    parseFastqFile(file, file_data, flip=flip)
    
    start_amino_barcode = aminoConversion_simplified( start_barcode)
    end_amino_barcode = aminoConversion_simplified(end_barcode)

    start_counts = 0
    end_counts = 0
    both_counts = 0

    start_counts_amino = 0
    end_counts_amino = 0
    both_counts_amino = 0

    for key, value in file_data.items():

        scores, seq = np.divmod(value["primarySequence"], 5)
            
        lowQualityCount = np.count_nonzero(np.less(scores, minQualityScore))
        
        if min(scores) == 0 :
            continue

        protein = aminoConversion(seq, **kwargs)
            
        stringRep = str(np.take(baseConversionArray, seq).tobytes())
        
        contain_start = start_barcode in stringRep
        contain_end = end_barcode in stringRep
        
        start_counts_amino += 1 if start_amino_barcode in protein else 0
        end_counts_amino += 1 if end_amino_barcode in protein else 0
        both_counts_amino += 1 if (start_amino_barcode in protein and end_amino_barcode in protein) else 0
        
        # print(f"Forward Read ID: {key}, Finalized Sequence: {converted_seq}")
        # print(f"Foward barcode CCCGGGTACCTTTCTATTCTCACTCTTCTTGT in finalized sequence: {'CCCGGGTACCTTTCTATTCTCACTCTTCTTGT' in converted_seq}")
        # print(f"End Barcode TGTGGTGGAGGTTCGGCCGGGCGCGGTGGT in finalized sequence: {'TGTGGTGGAGGTTCGGCCGGGCGCGGTGGT' in converted_seq}")
        
        start_counts += 1 if contain_start else 0
        end_counts += 1 if contain_end else 0
        both_counts += 1 if contain_start and contain_end else 0
       
    return {
        "total_reads": len(file_data),
        "start_counts": start_counts,
        "end_counts": end_counts,
        "both_counts": both_counts,
        "start_counts_amino": start_counts_amino,
        "end_counts_amino": end_counts_amino,
        "both_counts_amino": both_counts_amino
    }
    
def evalutate_json_file(file, flip = False, start_barcode = "AAA", end_barcode = "TTT", **kwargs):
    with open(file, 'r') as infile:
        file_data = json.load(infile)
    
    return evalutate_dict_file(file_data, flip=flip, start_barcode=start_barcode, end_barcode=end_barcode, **kwargs)
    
def evalutate_dict_file(file_data, flip = False, start_barcode = "AAA", end_barcode = "TTT", **kwargs):
    
    
    start_amino_barcode = aminoConversion_simplified( start_barcode)
    end_amino_barcode = aminoConversion_simplified(end_barcode)

    start_counts = 0
    end_counts = 0
    both_counts = 0

    start_counts_amino = 0
    end_counts_amino = 0
    both_counts_amino = 0

    for key, value in file_data.items():
        
        
        
        contain_start = start_barcode in value["sequences"]
        contain_end = end_barcode in value["sequences"]
        
        start_counts_amino += 1 if start_amino_barcode in value.get("proteinSequence", "") else 0
        end_counts_amino += 1 if end_amino_barcode in value.get("proteinSequence", "") else 0
        both_counts_amino += 1 if (start_amino_barcode in value.get("proteinSequence", "") and end_amino_barcode in value["proteinSequence"]) else 0
        
        # print(f"Forward Read ID: {key}, Finalized Sequence: {converted_seq}")
        # print(f"Foward barcode CCCGGGTACCTTTCTATTCTCACTCTTCTTGT in finalized sequence: {'CCCGGGTACCTTTCTATTCTCACTCTTCTTGT' in converted_seq}")
        # print(f"End Barcode TGTGGTGGAGGTTCGGCCGGGCGCGGTGGT in finalized sequence: {'TGTGGTGGAGGTTCGGCCGGGCGCGGTGGT' in converted_seq}")
        
        start_counts += 1 if contain_start else 0
        end_counts += 1 if contain_end else 0
        both_counts += 1 if contain_start and contain_end else 0
       
    return {
        "total_reads": len(file_data),
        "start_counts": start_counts,
        "end_counts": end_counts,
        "both_counts": both_counts,
        "start_counts_amino": start_counts_amino,
        "end_counts_amino": end_counts_amino,
        "both_counts_amino": both_counts_amino
    }
    


def _return_data_files():
    return [["data/GAL_CA_R11/PID-1309-GAL-CA-1-PC_S105_R1_001.fastq","data/GAL_CA_R11/PID-1309-GAL-CA-1-PC_S105_R2_001.fastq"],
            ["data/GAL_CA_R12/PID-1309-GAL-CA-2-PC_S106_R1_001.fastq","data/GAL_CA_R12/PID-1309-GAL-CA-2-PC_S106_R2_001.fastq"],
            ["data/GAL_CA_R13/PID-1309-GAL-CA3-PC_S91_R1_001.fastq","data/GAL_CA_R13/PID-1309-GAL-CA3-PC_S91_R2_001.fastq"],
            ["data/MON_BSA/PID-1309-M7-MON-BSA-1_S87_R1_001.fastq","data/MON_BSA/PID-1309-M7-MON-BSA-1_S87_R2_001.fastq"],
            ["data/MON_BSA/PID-1309-M7-MON-BSA-2_S88_R1_001.fastq","data/MON_BSA/PID-1309-M7-MON-BSA-1_S87_R2_001.fastq"],
            ["data/MON_CONA/PID-1309-M7-MON-CONA-1_S85_R1_001.fastq","data/MON_CONA/PID-1309-M7-MON-CONA-1_S85_R2_001.fastq"],
            ["data/MON_CONA/PID-1309-M7-MON-CONA-2_S86_R1_001.fastq", "data/MON_CONA/PID-1309-M7-MON-CONA-2_S86_R2_001.fastq"]
            ]
    
def evalutate_p3anut(start_barcode="TATTCTCACTCTTCT", end_barcode="GGTGGAGGTTCG", parse_args={}):
    files = _return_data_files()
    
    results = {}
    
    for forward, reverse in files:
        print(f"Evaluating files: {forward}, {reverse}")
        t1 = evalutate_fastq_file(forward, start_barcode=start_barcode, end_barcode=end_barcode, flip=False)
        t2 = evalutate_fastq_file(reverse, start_barcode=start_barcode, end_barcode=end_barcode, flip=True)
        
        merged_data = {}
        _ = parse(forward, reverse, data=merged_data, **parse_args)
        
        t3 = evalutate_dict_file(merged_data, flip=False, start_barcode=start_barcode, end_barcode=end_barcode)
        
        results[forward] = {
            "forward": t1,
            "reverse": t2,
            "merged": t3
        }
        
        print(f"Forward File Evaluation: {t1['both_counts']/t1['total_reads']*100:.8f}% reads contain both barcodes")
        print(f"Reverse File Evaluation: {t2['both_counts']/t2['total_reads']*100:.8f}% reads contain both barcodes")
        print(f"Merged File Evaluation: {t3['both_counts']/t3['total_reads']*100:.8f}% reads contain both barcodes")
        
    return results

def evalutate_delta(forfile, reversefile, d_start = 1, d_end = 19, start_barcode="GACTATTCTCACTCTTCT", end_barcode="GGTGGAGGTTCG",):
    results = {}
    
    for d_s in range(d_start, d_end + 1):
        parse_args = {"scoreOffset": d_s, "multiprocess": True, "cull_maxlength": 100}
        print(f"Evaluating with score offset: {d_s}")
        
        merged_data = {}
        _ = parse(forfile, reversefile, data=merged_data, **parse_args)
        
        results[d_s] = evalutate_dict_file(merged_data, flip=False, start_barcode=start_barcode, end_barcode=end_barcode)
        
            
    return results
    
if __name__ == "__main__":
    
    file_list = _return_data_files()    
    
    for forFile, revFile in file_list:
        print(os.path.exists(forFile), os.path.exists(revFile))
    
    t = evalutate_p3anut(start_barcode="TATTCTCACTCTTCT", end_barcode="GGTGGAGGTTCG", parse_args={"multiprocess": True, "cull_maxlength": 100, "scoreOffset":11})
    
    with open("dev_tools/p3anut_evaluation.json", 'w') as outfile:
        json.dump(t, outfile, indent=4)
        
        
    # t = evalutate_delta(file_list[0][0], file_list[0][1], d_start=1, d_end=19, start_barcode="GACTATTCTCACTCTTCT", end_barcode="GGTGGAGGTTCG")
    # with open("dev_tools/p3anut_delta_evaluation.json", 'w') as outfile:
    #     json.dump(t, outfile, indent=4)
    
    # forFile = "/Users/ethankoland/Desktop/Undergrad/Year 3/3rd Year Project/code/data/PID-1309-GAL-BSA-1-PC_S107_R1_001.fastq"
    # revFile = "/Users/ethankoland/Desktop/Undergrad/Year 3/3rd Year Project/code/data/PID-1309-GAL-BSA-1-PC_S107_R2_001.fastq"
    # mergeFile = "/Users/ethankoland/Desktop/Side Projects/P3ANUT/dev_Tools/mergeMismatchData.json"

    # t1 = evalutate_fastq_file(forFile, start_barcode="GACTATTCTCACTCTTCT", end_barcode="GGTGGAGGTTCG", flip=False, minQualityScore=20, minQualityScoreCount=5, proteinConversion=True, aminoBaseRange=6)
    # t2 = evalutate_fastq_file(revFile, start_barcode="GACTATTCTCACTCTTCT", end_barcode="GGTGGAGGTTCG", flip=True, minQualityScore=20, minQualityScoreCount=5, proteinConversion=True, aminoBaseRange=6)
    # t3 = evalutate_json_file(mergeFile, start_barcode="GACTATTCTCACTCTTCT", end_barcode="GGTGGAGGTTCG", flip=False, minQualityScore=20, minQualityScoreCount=5, proteinConversion=True, aminoBaseRange=6)

    # print(f"Forward File Evaluation: {t1["both_counts"]/t1["total_reads"]*100:.8f}% reads contain both barcodes")
    # print(f"Reverse File Evaluation: {t2["both_counts"]/t2["total_reads"]*100:.8f}% reads contain both barcodes")
    # print(f"Merged File Evaluation: {t3["both_counts"]/t3["total_reads"]*100:.8f}% reads contain both barcodes")
    
    
    
    # t = evalutate_delta(forFile, revFile, d_start=10, d_end=30)
    # json.dump(t, open("dev_tools/delta_evaluation.json", 'w'), indent=4)
    
    


# print(evalutate_fastq_file(forFile, start_barcode="GACTATTCTCACTCTTCTTGT", end_barcode="TGTGGTGGAGGTTCG", flip=False, minQualityScore=20, minQualityScoreCount=5, proteinConversion=True, aminoBaseRange=6))
# print(evalutate_fastq_file(revFile, start_barcode="GACTATTCTCACTCTTCTTGT", end_barcode="TGTGGTGGAGGTTCG", flip=True, minQualityScore=20, minQualityScoreCount=5, proteinConversion=True, aminoBaseRange=6))
# # print(evalutate_json_file(mergeFile, start_barcode="GACTATTCTCACTCTTCT", end_barcode="GGTGGAGGTTCG", flip=False, minQualityScore=20, minQualityScoreCount=5, proteinConversion=True, aminoBaseRange=6))
# # print(evalutate_json_file(mergeFile, start_barcode="GACTATTCTCACTCTTCT", end_barcode="TGTGGTGGAGGTTCG", flip=False, minQualityScore=20, minQualityScoreCount=5, proteinConversion=True, aminoBaseRange=6))
# print(evalutate_json_file(mergeFile, start_barcode="GACTATTCTCACTCTTCTTGT", end_barcode="TGTGGTGGAGGTTCG", flip=False, minQualityScore=20, minQualityScoreCount=5, proteinConversion=True, aminoBaseRange=6))

# counts = []
# data = []
# with open("data/PID-1309-GAL-CA-1-PC_S105_R1_001.fastq", 'r') as infile:
#     line_counter = 0
#     for line in infile:
#         line_counter += 1
#         if line_counter % 4 == 2:
#             data.append(line.strip())
#             counts.append(len(line.strip()))
            
# counts = np.array(counts)
# print(f"Mean: {np.mean(counts)}, Median: {np.median(counts)}, Std: {np.std(counts)}")
# unique, unique_counts = np.unique(counts, return_counts=True)
# for u, c in zip(unique, unique_counts):
#     print(f"Length: {u}, Count: {c}")