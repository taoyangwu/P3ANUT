import pandas as pd
import csv
import numpy as np
import json
import os
import re

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


def evalutate_fastq_file(file, flip = False, start_barcode = "AAA", end_barcode = "TTT", target_length = 68,
                         upsilion_statement = None, phi_statement = None, **kwargs):
    
    
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

    upsilon_score = 0
    phi_score = 0

    sequence_length_counts = {}

    for key, value in file_data.items():

        scores, seq = np.divmod(value["primarySequence"], 5)
            
        lowQualityCount = np.count_nonzero(np.less(scores, minQualityScore))

        if min(scores) == 0 :
            continue

        protein = aminoConversion(seq, **kwargs)
            
        stringRep = str(np.take(baseConversionArray, seq).tobytes())[2:-1]

        sequence_length_counts[len(stringRep)] = sequence_length_counts.get(len(stringRep), 0) + 1
        most_frequent_length = max(sequence_length_counts, key=sequence_length_counts.get)

        upsilon_score += 1 if upsilion_statement != None and re.match(upsilion_statement, stringRep) else 0
        phi_score += 1 if phi_statement != None and re.match(phi_statement, stringRep) else 0
    
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

    output = {
        "total_reads": len(file_data),
        "start_counts": start_counts,
        "end_counts": end_counts,
        "both_counts": both_counts,
        "start_counts_amino": start_counts_amino,
        "end_counts_amino": end_counts_amino,
        "both_counts_amino": both_counts_amino,
        "sequence_length_score" : sequence_length_counts.get(most_frequent_length, 0) / sum(sequence_length_counts.values())
    }

    if "upsilon_score" != None:
        output["upsilon_score"] = upsilon_score
    
    if "phi_score" != None:
        output["phi_score"] = phi_score

    return output





def calculate_overlap(file1, file2):
    def _load(file_path):
        sequences = set()

        cull_minlength = 8
        cull_maxlength = 256
        fileSeperator = r"\+"

        with open(file_path, "r") as file:
            #Example Entry after regex as a tuple
            #('NB501061:163:HVYLLAFX3:1:11101:1980:1063' - Run ID,  
            # '1:N:0:GTATTATCT+CATATCGTT' - Additional Run Information, 
            # 'TGTAGACTATTCTCACTCTTCTTGTCTGGTTCCTCCGCGTCCGACGTGTGGTGGAGGTTCGGTCGACG', - DNA Sequence 
            # 'AAAAAEE<EE<EEEEEEEEEEAEEEE/EEAEEEEA//AA/EAAEEEEEEEEAEEEEEE/EA<EEEEE6' - DNA Quality Score)
            regexExpression = r"@([A-Z0-9.:\-\s]+)(?:\s)([A-Z0-9:+\/-]+)(?:\s*)([CODONS]{cull_minlength,cull_maxlength})(?:\s+SPLIT\s*)([!-I]{cull_minlength,cull_maxlength})".replace("cull_minlength", str(cull_minlength)).replace("cull_maxlength", str(cull_maxlength))
            regexExpression = regexExpression.replace("CODONS", "ATGCN")
            
            cleanedFileSeparator = fileSeperator.replace("\\\\", "\\")
            regexExpression = regexExpression.replace("SPLIT", f"[{cleanedFileSeparator}]?")
            entries = re.findall(regexExpression, file.read())

        for entry in entries:
            sequence = entry[0]
            sequences.add(sequence)
        return sequences
    
    sequences1 = _load(file1)
    sequences2 = _load(file2)

    overlap = sequences1.union(sequences2)

    return len(overlap), len(overlap) - len(sequences1), len(overlap) - len(sequences2)

    
def evalutate_json_file(file, flip = False, start_barcode = "AAA", end_barcode = "TTT", **kwargs):
    with open(file, 'r') as infile:
        file_data = json.load(infile)
    
    return evalutate_dict_file(file_data, flip=flip, start_barcode=start_barcode, end_barcode=end_barcode, **kwargs)
    
def evalutate_dict_file(file_data, flip = False, start_barcode = "AAA", end_barcode = "TTT", target_sequence_length = 68,
                        upsilion_statement = None, phi_statement = None, **kwargs):
    
    
    start_amino_barcode = aminoConversion_simplified( start_barcode)
    end_amino_barcode = aminoConversion_simplified(end_barcode)

    start_counts = 0
    end_counts = 0
    both_counts = 0

    start_counts_amino = 0
    end_counts_amino = 0
    both_counts_amino = 0

    upsilon_score = 0
    phi_score = 0

    sequences_lengths = {}

    for key, value in file_data.items():
        
        sequence_length = len(value["sequences"])
        
        contain_start = start_barcode in value["sequences"]
        contain_end = end_barcode in value["sequences"]
        
        start_counts_amino += 1 if start_amino_barcode in value.get("proteinSequence", "") else 0
        end_counts_amino += 1 if end_amino_barcode in value.get("proteinSequence", "") else 0
        both_counts_amino += 1 if (start_amino_barcode in value.get("proteinSequence", "") and end_amino_barcode in value["proteinSequence"]) else 0
        
        sequences_lengths[sequence_length] = sequences_lengths.get(sequence_length, 0) + 1

        upsilon_score += 1 if upsilion_statement != None and re.match(upsilion_statement, value["sequences"]) else 0
        phi_score += 1 if phi_statement != None and re.match(phi_statement, value["sequences"]) else 0

        # print(f"Forward Read ID: {key}, Finalized Sequence: {converted_seq}")
        # print(f"Foward barcode CCCGGGTACCTTTCTATTCTCACTCTTCTTGT in finalized sequence: {'CCCGGGTACCTTTCTATTCTCACTCTTCTTGT' in converted_seq}")
        # print(f"End Barcode TGTGGTGGAGGTTCGGCCGGGCGCGGTGGT in finalized sequence: {'TGTGGTGGAGGTTCGGCCGGGCGCGGTGGT' in converted_seq}")
        
        start_counts += 1 if contain_start else 0
        end_counts += 1 if contain_end else 0
        both_counts += 1 if contain_start and contain_end else 0
       
    output = {
        "total_reads": len(file_data),
        "start_counts": start_counts,
        "end_counts": end_counts,
        "both_counts": both_counts,
        "start_counts_amino": start_counts_amino,
        "end_counts_amino": end_counts_amino,
        "both_counts_amino": both_counts_amino,
        "sequence_length_score" : sequences_lengths.get(target_sequence_length, 0) / sum(sequences_lengths.values())
    }

    if "upsilon_score" != None:
        output["upsilon_score"] = upsilon_score
    
    if "phi_score" != None:
        output["phi_score"] = phi_score

    return output
    


def _return_data_files():
    return [["data/GAL_CA_R11/PID-1309-GAL-CA-1-PC_S105_R1_001.fastq","data/GAL_CA_R11/PID-1309-GAL-CA-1-PC_S105_R2_001.fastq"],
            ["data/GAL_CA_R12/PID-1309-GAL-CA-2-PC_S106_R1_001.fastq","data/GAL_CA_R12/PID-1309-GAL-CA-2-PC_S106_R2_001.fastq"],
            ["data/GAL_CA_R13/PID-1309-GAL-CA3-PC_S91_R1_001.fastq","data/GAL_CA_R13/PID-1309-GAL-CA3-PC_S91_R2_001.fastq"],
            ["data/MON_BSA/PID-1309-M7-MON-BSA-1_S87_R1_001.fastq","data/MON_BSA/PID-1309-M7-MON-BSA-1_S87_R2_001.fastq"],
            ["data/MON_BSA/PID-1309-M7-MON-BSA-2_S88_R1_001.fastq","data/MON_BSA/PID-1309-M7-MON-BSA-1_S87_R2_001.fastq"],
            ["data/MON_CONA/PID-1309-M7-MON-CONA-1_S85_R1_001.fastq","data/MON_CONA/PID-1309-M7-MON-CONA-1_S85_R2_001.fastq"],
            ["data/MON_CONA/PID-1309-M7-MON-CONA-2_S86_R1_001.fastq", "data/MON_CONA/PID-1309-M7-MON-CONA-2_S86_R2_001.fastq"]
            ]

def _large_evaluation(file_list):
    pairs = []
    with open(file_list, 'r') as f:
        lines = f.readlines()
        pairs = [line.strip().split() for line in lines]

    return pairs
    
    
def evalutate_p3anut(files, start_barcode="TATTCTCACTCTTCT", end_barcode="GGTGGAGGTTCG", parse_args={}):
    # files = _return_data_files()
    
    results = {}

    upsilion_statement = "[TAGC]{7}start_barcode[TGCA]{27}end_barcode[TAGC]{7}".replace("start_barcode", start_barcode).replace("end_barcode", end_barcode)
    phi_statement = "[TAGC]*start_barcode[TGCA]{27}end_barcode[TAGC]*".replace("start_barcode", start_barcode).replace("end_barcode", end_barcode)
    
    for forward, reverse in files:
        print(f"Evaluating files: {forward}, {reverse}")

        t1 = evalutate_fastq_file(forward, start_barcode=start_barcode, end_barcode=end_barcode, flip=False, upsilion_statement=upsilion_statement, phi_statement=phi_statement, **parse_args)

        try:
            t1 = evalutate_fastq_file(forward, start_barcode=start_barcode, end_barcode=end_barcode, flip=False, upsilion_statement=upsilion_statement, phi_statement=phi_statement, **parse_args)
            t2 = evalutate_fastq_file(reverse, start_barcode=start_barcode, end_barcode=end_barcode, flip=True, upsilion_statement=upsilion_statement, phi_statement=phi_statement, **parse_args)
        except Exception as e:
            with open("error_log.txt", 'a') as error_file:
                error_file.write(f"Error evaluating files {forward} and {reverse}: {e}\n")
            continue

        sequence_length_count, f1_count, f2_count = calculate_overlap(forward, reverse)
        
        merged_data = {}
        t0 = parse(forward, reverse, data=merged_data, **parse_args)
        
        t1["retention_rate"] = (sequence_length_count - f1_count) / sequence_length_count
        t2["retention_rate"] = (sequence_length_count - f2_count) / sequence_length_count

        t3 = evalutate_dict_file(merged_data, flip=False, start_barcode=start_barcode, end_barcode=end_barcode, upsilion_statement=upsilion_statement, phi_statement=phi_statement, **parse_args)

        t3["retention_rate"] = t3["total_reads"] / sequence_length_count if sequence_length_count > 0 else 0
        
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

def evaluate_folder(path, start_barcode="GACTATTCTCACTCTTCT", end_barcode="GGTGGAGGTTCG", **kwargs):
    results = {}
    
    for file in os.listdir(path):
        if file.endswith(".fastq"):
            file_path = os.path.join(path, file)
            print(f"Evaluating file: {file_path}")
            results[file] = evalutate_fastq_file(file_path, start_barcode="GACTATTCTCACTCTTCT", end_barcode="GGTGGAGGTTCG", flip=False)
    
    return results

def calculate_file_lengths(file_list):
    lengths = {}

    run_name_lambda = lambda file_path: os.path.basename(file_path).split("_R1_")[0]

    def _load(file_path):
        sequences = set()

        cull_minlength = 8
        cull_maxlength = 256
        fileSeperator = r"\+"

        with open(file_path, "r") as file:
            #Example Entry after regex as a tuple
            #('NB501061:163:HVYLLAFX3:1:11101:1980:1063' - Run ID,  
            # '1:N:0:GTATTATCT+CATATCGTT' - Additional Run Information, 
            # 'TGTAGACTATTCTCACTCTTCTTGTCTGGTTCCTCCGCGTCCGACGTGTGGTGGAGGTTCGGTCGACG', - DNA Sequence 
            # 'AAAAAEE<EE<EEEEEEEEEEAEEEE/EEAEEEEA//AA/EAAEEEEEEEEAEEEEEE/EA<EEEEE6' - DNA Quality Score)
            regexExpression = r"@([A-Z0-9.:\-\s]+)(?:\s)([A-Z0-9:+\/-]+)(?:\s*)([CODONS]{cull_minlength,cull_maxlength})(?:\s+SPLIT\s*)([!-I]{cull_minlength,cull_maxlength})".replace("cull_minlength", str(cull_minlength)).replace("cull_maxlength", str(cull_maxlength))
            regexExpression = regexExpression.replace("CODONS", "ATGCN")
            
            cleanedFileSeparator = fileSeperator.replace("\\\\", "\\")
            regexExpression = regexExpression.replace("SPLIT", f"[{cleanedFileSeparator}]?")
            entries = re.findall(regexExpression, file.read())

        return len(entries)
    
    for forward, reverse in file_list:

        run_name = run_name_lambda(forward)

        try:
            f_len = _load(forward)
            r_len = _load(reverse)
        except Exception as e:
            print(f"Error calculating lengths for {forward} and {reverse}: {e}")
            continue

        lengths[run_name] = {
            "forward_length": f_len,
            "reverse_length": r_len,
            "total_length": f_len + r_len
        }

    return lengths

def t():
    file_list = [["/home/proxima/Desktop/Side_Projects/P3ANUT/data/Original_Sequencing/M7-Gal/PID-1309-GAL-CA3-PC_S91_R1_001.fastq", 
                 "/home/proxima/Desktop/Side_Projects/P3ANUT/data/Original_Sequencing/M7-Gal/PID-1309-GAL-CA3-PC_S91_R2_001.fastq"]]
    
    for forFile, revFile in file_list:
        print(os.path.exists(forFile), os.path.exists(revFile))
    
    t = evalutate_p3anut(file_list, start_barcode="TATTCTCACTCTTCT", end_barcode="GGTGGAGGTTCG", parse_args={"multiprocess": True, "cull_maxlength": 100, "scoreOffset":11})
    
    with open("p3anut_evaluation_ups_phi_2.json", 'w') as outfile:
        json.dump(t, outfile, indent=4)


    
if __name__ == "__main__":
    
    # file_list = _return_data_files() 
    # 
    t()   

    file_list = _large_evaluation("/home/proxima/Desktop/Side_Projects/FLASH_CASPAR/missed_files.txt")

    file_lengths = calculate_file_lengths(file_list)

    with open("file_lengths_2.json", 'w') as outfile:
        json.dump(file_lengths, outfile, indent=4)
    
    # t = evalutate_p3anut(file_list, start_barcode="TATTCTCACTCTTCT", end_barcode="GGTGGAGGTTCG", parse_args={"multiprocess": True, "cull_maxlength": 100, "scoreOffset":11})
    
    # with open("p3anut_evaluation_ups_phi.json", 'w') as outfile:
    #     json.dump(t, outfile, indent=4)
        
        
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