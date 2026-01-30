import os
from math import tanh, floor
import json
from shutil import rmtree
import matplotlib.pyplot as plt
import numpy as np
from math import inf
import re

def split_dataset(filePath, barcodes, save_Path, fileBase, start_code="ATATTTATG",
                   end_code="TGCGGTGGA", key_barcode_length=9, hamming_correction=True):
    '''
    @LH00157:478:22MH22LT4:6:1101:9683:1098 1:N:0:ANTGAGCG+TGGAGAGA
    ATACATACTAATACGACTCACTATAGGATTAAGGAGGTGATATTTATGCATCCGGGCCATCTGAAAGGCGCGAAATTGGCATGTGGTATGCGAAAAAACAGGTTAATGTTAAGTGTTGTTCTTCATGGTAATTTTTATTGCGGTGGAGGA
    +
    IIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII99II9II-9IIIIIIIII9IIIIIII9II9IIIII-IIIIIIIIIIIIIIIIIII9I9III9IIIIIIII-99IIIIII9IIIII-99II9-II
    '''
    
    create_save_path = lambda path, filebase, prefix: os.path.join(path, f"{prefix}_{filebase}")
    # extact_bases = {}
    
    # reference_dict = {}
    # file_Counts = {}
    
    # for key, value in barcodes.items():
        
    #     local_path = create_save_path(save_Path, fileBase, key)
    #     os.makedirs(local_path, exist_ok=True)
    #     reference_dict[value] = local_path
    #     file_Counts[key] = 0
        
        
    
    # os.makedirs(local_path, exist_ok=True)
    # non_matched = local_path
    # file_Counts["non_matched"] = 0

    reference_dict, file_Counts = _make_reference_dicts(barcodes, save_Path, fileBase)
    non_matched_path = create_save_path(save_Path, fileBase, "non_matched")

    start_positions, end_positions, total_lenghts = [], [], []
    temp_score = 0
        
    with open(filePath, 'r') as in_file:
        buffer = []
        for line in in_file:
            buffer.append(line.strip())
            if len(buffer) == 4:
                header, sequence, plus, quality = buffer
                barcode_seq = sequence[0:4].upper()
                start_cut, end_cut = _cut_sequence_front_biasis(sequence.upper(), start_code, end_code)
                
                
                start_positions.append(start_cut)
                end_positions.append(end_cut)
                total_lenghts.append(len(sequence))

                if barcode_seq in reference_dict:
                    save_dir = reference_dict[barcode_seq]
                    file_Counts[barcode_seq] = file_Counts.get(barcode_seq, 0) + 1
                    with open(os.path.join(save_dir, f"{fileBase}.fastq"), 'a') as out_file:
                        out_file.write(f"{header}\n{sequence[start_cut : end_cut]}\n{plus}\n{quality[start_cut : end_cut]}\n")
                else:

                    barcode_seq_corrected = (_hamming_corrections(barcode_seq, barcodes, max_distance=1)
                                                if hamming_correction else "non_matched" )
                    temp_score += 1 if barcode_seq_corrected != "non_matched" else 0
                    save_dir = reference_dict[barcode_seq_corrected] 
                    file_Counts[barcode_seq_corrected] = file_Counts.get(barcode_seq_corrected, 0) + 1
                    with open(os.path.join(save_dir, f"{fileBase}.fastq"), 'a') as out_file:
                        out_file.write(f"{header}\n{sequence[start_cut : end_cut]}\n{plus}\n{quality[start_cut : end_cut]}\n")
                
                buffer = []  # Clear buffer for next record
                
    print("File split complete. Summary:")
    print(f"Hamming corrections made: {temp_score}")
    for key, count in file_Counts.items():
        print(f"{key}: {count} reads")
    
    return file_Counts, [start_positions, end_positions, total_lenghts]

def _make_reference_dicts(barcodes, save_Path, fileBase):
    reference_dict = {}
    file_counts = {}

    modified_barcodes = barcodes | {"non_matched": "non_matched"}

    for key, value in modified_barcodes.items():
        local_path = os.path.join(save_Path, f"{key}_{fileBase}")
        os.makedirs(local_path, exist_ok=True)
        reference_dict[value] = local_path
        file_counts[key] = 0
    return reference_dict, file_counts

def _hamming_corrections(barcode_seq, barcodes, max_distance=1):
    values = {}
    min_distance = inf

    for key, value in barcodes.items():
        distance = sum(c1 != c2 for c1, c2 in zip(barcode_seq, value))
        values.setdefault(distance, []).append(value)
        min_distance = min(min_distance, distance)

    if min_distance <= max_distance and len(values[min_distance]) == 1:
        return values[min_distance][0]
    else:
        return "non_matched"

def two_file_split_dataset(filePath1, filePath2, f1_barcodes, f2_barcodes, 
                           save_Path, fileBase, start_code="ATATTTATG",
                   end_code="TGCGGTGGA", key_barcode_length=9, ):
    
    pass

def analyze_file_bardcode_overlap(file1Path, file2Path, f1_barcodes, f2_barcodes):
    # barcodes are "R3" : "CCAT", ...
    f1_inverted = {value: key for key, value in f1_barcodes.items()} | {"non_matched": "non_matched"}
    f2_inverted = {value: key for key, value in f2_barcodes.items()} | {"non_matched": "non_matched"}
    


    file1_sets = {}
    file1_qualityScores = np.zeros((4, 41))  # 4 bases, quality scores 0-40

    phred_lambda = lambda q: ord(q) - 33

    with open(file1Path, 'r') as in_file:
        buffer = []
        for line in in_file:
            buffer.append(line.strip())
            if len(buffer) == 4:
                header, sequence, plus, quality = buffer
                run_id = header.split(' ')[0]
                barcode_seq = sequence[0:4].upper()

                if barcode_seq in f1_inverted:
                    barcode_key = f1_inverted[barcode_seq]
                else:
                    barcode_key = "non_matched"
                
                file1_sets[run_id] = barcode_key

                for i in range(4):
                    q_score = phred_lambda(quality[i])
                    file1_qualityScores[i][q_score] += 1

                buffer = []  # Clear buffer for next record

    file2_sets = {}
    file2_qualityScores = np.zeros((4, 41))  # 4 bases, quality scores 0-40

    with open(file2Path, 'r') as in_file:
        buffer = []
        for line in in_file:
            buffer.append(line.strip())
            if len(buffer) == 4:
                header, sequence, plus, quality = buffer
                run_id = header.split(' ')[0]
                barcode_seq = sequence[0:4].upper()

                if barcode_seq in f2_inverted:
                    barcode_key = f2_inverted[barcode_seq]
                else:
                    barcode_key = "non_matched"
                
                file2_sets[run_id] = barcode_key

                for i in range(4):
                    q_score = phred_lambda(quality[i])
                    file2_qualityScores[i][q_score] += 1

                buffer = []  # Clear buffer for next record

    f1_ref_index = f1_barcodes| {"non_matched"  : "non_matched", "not_found" : "not_found"}
    f1_index = {}
    for i, key in enumerate(f1_ref_index.keys()):
        f1_index[key] = i

    f2_ref_index = f2_barcodes| {"non_matched"  : "non_matched", "not_found" : "not_found"}
    f2_index = {}
    for i, key in enumerate(f2_ref_index.keys()):
        f2_index[key] = i

    overlap_matrix = np.zeros((len(f1_index), len(f2_index)), dtype=int)

    run_ids =set(file1_sets.keys()).union(set(file2_sets.keys()))

    for run_id in run_ids:
        f1_barcode = file1_sets.get(run_id, "not_found")
        f2_barcode = file2_sets.get(run_id, "not_found")

        f1_idx = f1_index[f1_barcode]
        f2_idx = f2_index[f2_barcode]

        overlap_matrix[f1_idx][f2_idx] += 1

    




    diagonal_sum = np.trace(overlap_matrix)
    total_counts = np.sum(overlap_matrix)
    print(f"Diagonal Sum (Matching Barcodes): {diagonal_sum}")
    print(f"Total Counts: {total_counts}")
    print(f"Percentage of Matching Barcodes: {diagonal_sum / total_counts * 100:.2f}%")

    f1_qs_nonZero = np.nonzero(file1_qualityScores)
    f2_qs_nonZero = np.nonzero(file2_qualityScores)

    for i in range(len(f1_qs_nonZero[0])):
        base_idx = f1_qs_nonZero[0][i]
        score_idx = f1_qs_nonZero[1][i]
        count = file1_qualityScores[base_idx][score_idx]
        print(f"File 1 - Base {base_idx + 1}, Quality Score {score_idx}: {int(count)} reads")

    #2 columns
    figure, [ax, ax1] = plt.subplots(ncols=2, figsize=(16, 8))

    
    for (i, j), val in np.ndenumerate(overlap_matrix):
        ax.text(j, i, val, ha='center', va='center', color='black')
    cax = ax.matshow(overlap_matrix, cmap='Blues')
    plt.colorbar(cax)
    ax.set_xticklabels([''] + list(f2_index.keys()))
    ax.set_yticklabels([''] + list(f1_index.keys()))
    ax.set_xlabel('File 2 Barcodes')
    ax.set_ylabel('File 1 Barcodes')
    ax.set_title('Barcode Overlap Matrix')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis='x', rotation=45)

    #two bar plots for quality scores
    x = np.arange(41)
    width = 0.35  # the width of the bars
    bottom_f1 = np.zeros(41)
    bottom_f2 = np.zeros(41)
    for i in range(4):
        ax1.bar(x - width/2, file1_qualityScores[i], width, label=f'File 1 - Base {i+1}', bottom = bottom_f1)
        ax1.bar(x + width/2, file2_qualityScores[i], width, label=f'File 2 - Base {i+1}', bottom = bottom_f2)
        bottom_f1 += file1_qualityScores[i]
        bottom_f2 += file2_qualityScores[i]

    ax1.set_xlabel('Quality Score')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Quality Score Distribution for First 4 Bases')
    ax1.legend()

    plt.tight_layout()
    plt.show()






    plt.savefig("data/Rhau/barcode_overlap_matrix.png")


                
            






def _cut_sequence(sequence, start_barcode, end_barcode, barcode_length=9):
    
    max_barcode_length = max(len(start_barcode), len(end_barcode))
    
    padded_seq = "N" * max_barcode_length + sequence + "N" * max_barcode_length
    
    score_offset_multiplier = 4
    score_offset = [floor(tanh((i - len(padded_seq)/2) * (score_offset_multiplier / len(padded_seq))) * max_barcode_length) for i in range(len(padded_seq))]
    
    matched_bases_start = []
    matched_bases_start_end = []
    
    for i in range(len(sequence) + max_barcode_length):
        start_count = 0
        for j in range(len(start_barcode)):
            start_count += (1 if padded_seq[i + j] == start_barcode[j] else 0)
        matched_bases_start.append(start_count)
        end_count = 0
        for j in range(len(end_barcode)):
            end_count += (1 if padded_seq[i + j] == end_barcode[j] else 0)
        matched_bases_start_end.append(end_count)
    
    
    
    start_index = matched_bases_start.index(max(matched_bases_start))
    end_index = matched_bases_start_end[::-1].index(max(matched_bases_start_end))
    
    # t1 = [matched_bases_start[i] - score_offset[i] for i in range(len(matched_bases_start))]
    # t2 = [matched_bases_start_end[i] + score_offset[i] for i in range(len(matched_bases_start_end))]
    
    # start_index_2 = t1.index(max(t1)) - max_barcode_length
    # end_index_2 = t2[::-1].index(max(t2)) - max_barcode_length

    start_index_offset = start_index - barcode_length -1
    end_index_offset = end_index + barcode_length
    
    return start_index_offset,  end_index_offset

def _cut_sequence_front_biasis(sequence, start_barcode, end_barcode, length_of_random_sequence = 90,
                               barcode_length=9):
    
    matched_bases_start = []
   
    for i in range(len(sequence) - length_of_random_sequence - barcode_length):
        start_count = 0
        for j in range(len(start_barcode)):
            start_count += (1 if sequence[i + j] == start_barcode[j] else 0)
        matched_bases_start.append(start_count)

    start_index = matched_bases_start.index(max(matched_bases_start)) + len(start_barcode) - barcode_length

    end_center = start_index + length_of_random_sequence + barcode_length
    if end_center + barcode_length > len(sequence):
        return start_index + (len(start_barcode) - barcode_length), len(sequence)

    modified_end_barcode = end_barcode[:barcode_length]
    matched_end_bases = []

    for i in range(end_center - barcode_length, end_center + barcode_length):
        end_count = 0
        for j in range(min(barcode_length, len(sequence) - i)):
            end_count += (1 if sequence[i + j] == modified_end_barcode[j] else 0)
        matched_end_bases.append(end_count)


    end_index = matched_end_bases.index(max(matched_end_bases)) + end_center

    t1 = sequence[:start_index]
    t2 = sequence[end_index:]
    t_middle = sequence[start_index:end_index]
    t_middle_len = len(t_middle)


    return start_index, end_index

def graph_barcode_counts(counts_dict, title, save_path=None):
    
    
    labels = list(counts_dict.keys())
    counts = list(counts_dict.values())
    
    plt.figure(figsize=(10, 6))
    plt.bar(labels, counts, color='skyblue')
    plt.xlabel('Barcodes')
    plt.ylabel('Read Counts')
    plt.title(title)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(save_path if save_path else f"{title.replace(' ', '_')}.png")

def graph_total_length_distribution(lengths, save_path=None):

    max_length = max(lengths[2]) + 1
    
    start_bins = np.zeros((max_length))
    end_bins = np.zeros((max_length))
    random_bins = np.zeros((max_length))

    for start, end, total in zip(lengths[0], lengths[1], lengths[2]):
        start_bins[:start] += 1
        end_bins[end : total] += 1
        random_bins[start : end] += 1
    
    y = np.vstack([start_bins, random_bins, end_bins])
    labels = ['Start Segment', 'Random Segment', 'End Segment']
    colors = ['lightblue', 'lightgreen', 'salmon']

    plt.figure(figsize=(12, 5))
    plt.stackplot(range(max_length), y, labels=labels, colors=colors)
    plt.xlabel('Sequence Length')
    plt.ylabel('Frequency')
    plt.title('Total Length Distribution')
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(save_path if save_path else "total_length_distribution.png")
    
       
def graph_cut_positions(start_positions, end_positions, save_path=None):


    
    plt.figure(figsize=(12, 5))
    
    plt.subplot(1, 2, 1)
    plt.hist(start_positions, bins=50, color='lightgreen', edgecolor='black')
    plt.title('Start Cut Positions Distribution')
    plt.xlabel('Position')
    plt.ylabel('Frequency')
    
    plt.subplot(1, 2, 2)
    plt.hist(end_positions, bins=50, color='salmon', edgecolor='black')
    plt.title('End Cut Positions Distribution')
    plt.xlabel('Position')
    plt.ylabel('Frequency')
    
    plt.tight_layout()
    plt.savefig(save_path if save_path else "cut_positions_distribution.png")   




if __name__ == "__main__":
    '''
    R3 : ccat
    R4 : ccaa
    R5 : cgat
    R6 : cgtt
    R7 : atac
    '''

    invert_seq = lambda seq: ''.join({'A':'T', 'T':'A', 'C':'G', 'G':'C'}.get(base, base) for base in reversed(seq))
    
    barcode_forward = {
        'R3' : 'CCAT',
        'R4' : 'CCAA',
        'R5' : 'CGAT',
        'R6' : 'CGTT',
        'R7' : 'ATAC'
    }
    start_code = "ATACTAATACGACTCACTATAGGATTAAGGAGGTGATATTTATG"
    end_code = "TGCGGTGGAGGAGGAGGTAGCTAGGGACGGGGGGCGGGAGGCGGG" 

    # file_path = "data/Rhau/Raw_Data/Rhau18_12aa-R3-R7_R1.fq"
    # save_path = "data/Rhau/Forward_3"
    # rmtree(save_path, ignore_errors=True)
    

    # t1 = invert_seq(start_code)
    # t1_len = len(t1)
    # t2 = invert_seq(end_code)
    # t2_len = len(t2)

    


    # forward_counts, positions = split_dataset(file_path, barcode_forward, save_path, "Rhau18_12aa_F",
    #                                                                     start_code=start_code, end_code=end_code)
    # forward_counts = {"R3" : 1387953, "R4" : 1156818, "R5" : 1247606, "R6" : 866691, "R7" : 961845, "non_matched" : 163993}

    # print(max(positions[0]), min(positions[0]))
    # print(max(positions[1]), min(positions[1]))

    # forward_start = {}
    # for i in positions[0]:
    #     forward_start[i] = forward_start.get(i, 0) + 1

    # json.dump(forward_start, open("data/Rhau/Forward_2/start_positions_forward.json", 'w'), indent=4)

    # end_start = {}
    # for i in positions[1]:
    #     end_start[i] = end_start.get(i, 0) + 1 

    # json.dump(end_start, open("data/Rhau/Forward_2/end_positions_forward.json", 'w'), indent=4)

    # graph_total_length_distribution(positions, "data/Rhau/Forward_2/total_length_distribution_forward.png")
    # graph_barcode_counts(forward_counts, "Forward Read Barcode Distribution", "data/Rhau/Forward_2/barcode_distribution_forward.png")
    # graph_cut_positions(positions[0], positions[1], "data/Rhau/Forward_2/cut_positions_distribution_forward.png")
    
    '''
    R3 : ATCC (tagg)
    R4 : TTCC (aagg)
    R5 : ATGC (tacg)
    R6 : AAGC (ttcg)
    R7 : CTTA (gaat)
    R8 : GTTA (caat)
    '''
    # barcode_reverse = {
    #     "R3" : "CCTA",
    #     "R4" : "CCTT",
    #     "R5" : "CGTA",
    #     "R6" : "CGAA",
    #     "R7" : "ATTC",
    #     "R8" : "ATTG"
    # }

    barcode_reverse = {
        "R3" : "CCTA",
        "R4" : "CCTT",
        "R5" : "CGTA",
        "R6" : "CGAA",
        "R7" : "ATTC",
    }

    analyze_file_bardcode_overlap(
        "data/Rhau/Raw_Data/Rhau18_12aa-R3-R7_R1.fq",
        "data/Rhau/Raw_Data/Rhau18_12aa-R3-R7_R2.fq",
        barcode_forward,
        barcode_reverse)
    
    exit
    file_path = "data/Rhau/Raw_Data/Rhau18_12aa-R3-R7_R2.fq"
    save_path = "data/Rhau/Reverse"
    rmtree(save_path, ignore_errors=True)
    
    reverse_counts, positions = split_dataset(file_path, barcode_reverse, save_path, "Rhau18_12aa_R",
                                   invert_seq(end_code), invert_seq(start_code))
    
    forward_start = {}
    for i in positions[0]:
        forward_start[i] = forward_start.get(i, 0) + 1

    json.dump(forward_start, open("data/Rhau/Reverse/start_positions_reverse.json", 'w'), indent=4)

    end_start = {}
    for i in positions[1]:
        end_start[i] = end_start.get(i, 0) + 1 
    json.dump(end_start, open("data/Rhau/Reverse/end_positions_reverse.json", 'w'), indent=4)


    reverse_counts = {"R3" : 1407261, "R4" : 1020586, "R5" : 1220302, "R6" : 948704, "R7" : 957320, "R8" : 110, "non_matched" : 230623}
    graph_total_length_distribution(positions, "data/Rhau/Reverse/total_length_distribution_reverse.png")
    graph_barcode_counts(reverse_counts, "Reverse Read Barcode Distribution", "data/Rhau/Reverse/barcode_distribution_reverse.png")
    graph_cut_positions(positions[0], positions[1], "data/Rhau/Reverse/cut_positions_distribution_reverse.png")





