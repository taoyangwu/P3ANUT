import os
from math import tanh, floor
import json

def split_dataset(filePath, barcodes, save_Path, fileBase, start_code="ATATTTATG", end_code="TGCGGTGGA", key_barcode_length=9):
    '''
    @LH00157:478:22MH22LT4:6:1101:9683:1098 1:N:0:ANTGAGCG+TGGAGAGA
    ATACATACTAATACGACTCACTATAGGATTAAGGAGGTGATATTTATGCATCCGGGCCATCTGAAAGGCGCGAAATTGGCATGTGGTATGCGAAAAAACAGGTTAATGTTAAGTGTTGTTCTTCATGGTAATTTTTATTGCGGTGGAGGA
    +
    IIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII99II9II-9IIIIIIIII9IIIIIII9II9IIIII-IIIIIIIIIIIIIIIIIII9I9III9IIIIIIII-99IIIIII9IIIII-99II9-II
    '''
    
    create_save_path = lambda path, filebase, prefix: os.path.join(path, f"{prefix}_{filebase}")
    extact_bases = {}
    
    reference_dict = {}
    file_Counts = {}
    
    for key, value in barcodes.items():
        
        local_path = create_save_path(save_Path, fileBase, key)
        os.makedirs(local_path, exist_ok=True)
        reference_dict[value] = local_path
        file_Counts[key] = 0
        
        
    local_path = create_save_path(save_Path, fileBase, "non_matched")
    os.makedirs(local_path, exist_ok=True)
    non_matched = local_path
    file_Counts["non_matched"] = 0

    start_positions, end_positions = [], []
    
        
    with open(filePath, 'r') as in_file:
        buffer = []
        for line in in_file:
            buffer.append(line.strip())
            if len(buffer) == 4:
                header, sequence, plus, quality = buffer
                barcode_seq = sequence[0:4].upper()
                start_cut, end_cut = _cut_sequence(sequence.upper(), start_code, end_code)
                
                start_positions.append(start_cut)
                end_positions.append(end_cut)

                if barcode_seq in reference_dict:
                    save_dir = reference_dict[barcode_seq]
                    file_Counts[barcode_seq] = file_Counts.get(barcode_seq, 0) + 1
                    with open(os.path.join(save_dir, f"{fileBase}.fastq"), 'a') as out_file:
                        out_file.write(f"{header}\n{sequence[start_cut : end_cut]}\n{plus}\n{quality[start_cut : end_cut]}\n")
                else:
                    file_Counts["non_matched"] = file_Counts.get("non_matched", 0) + 1
                    with open(os.path.join(non_matched, f"{fileBase}.fastq"), 'a') as out_file:
                        out_file.write(f"{header}\n{sequence[start_cut : end_cut]}\n{plus}\n{quality[start_cut : end_cut]}\n")
                
                buffer = []  # Clear buffer for next record
                
    print("File split complete. Summary:")
    for key, count in file_Counts.items():
        print(f"{key}: {count} reads")
    
    return file_Counts, start_positions, end_positions

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

def graph_barcode_counts(counts_dict, title, save_path=None):
    import matplotlib.pyplot as plt
    
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
        
def graph_cut_positions(start_positions, end_positions, save_path=None):
    import matplotlib.pyplot as plt


    
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
    file_path = "data/Rhau/Raw_Data/Rhau18_12aa-R3-R7_R1.fq"
    save_path = "data/Rhau/Forward_2"
    start_code = "ATACTAATACGACTCACTATAGGATTAAGGAGGTGATATTTATG"
    end_code = "TGCGGTGGAGGAGGAGGTAGCTAGGGACGGGGGGCGGGAGGCGGG" 

    t1 = invert_seq(start_code)
    t2 = invert_seq(end_code)


    forward_counts, start_positions, end_positions = split_dataset(file_path, barcode_forward, save_path, "Rhau18_12aa_F",
                                                                        start_code=start_code, end_code=end_code)
    forward_counts = {"R3" : 1387953, "R4" : 1156818, "R5" : 1247606, "R6" : 866691, "R7" : 961845, "non_matched" : 163993}

    forward_start = {}
    for i in start_positions:
        forward_start[i] = forward_start.get(i, 0) + 1

    json.dump(forward_start, open("data/Rhau/Forward_2/start_positions_forward.json", 'w'), indent=4)

    end_start = {}
    for i in end_positions:
        end_start[i] = end_start.get(i, 0) + 1 

    json.dump(end_start, open("data/Rhau/Forward_2/end_positions_forward.json", 'w'), indent=4)

    graph_barcode_counts(forward_counts, "Forward Read Barcode Distribution", "data/Rhau/Forward_2/barcode_distribution_forward.png")
    graph_cut_positions(start_positions, end_positions, "data/Rhau/Forward_2/cut_positions_distribution_forward.png")
    
    '''
    R3 : ATCC (tagg)
    R4 : TTCC (aagg)
    R5 : ATGC (tacg)
    R6 : AAGC (ttcg)
    R7 : CTTA (gaat)
    R8 : GTTA (caat)
    '''
    barcode_reverse = {
        "R3" : "CCTA",
        "R4" : "CCTT",
        "R5" : "CGTA",
        "R6" : "CGAA",
        "R7" : "ATTC",
        "R8" : "ATTG"
    }
    
    file_path = "data/Rhau/Raw_Data/Rhau18_12aa-R3-R7_R2.fq"
    save_path = "data/Rhau/Reverse"
    
    reverse_counts, start_positions, end_positions = split_dataset(file_path, barcode_reverse, save_path, "Rhau18_12aa_R",
                                   invert_seq(start_code), invert_seq(end_code))
    
    forward_start = {}
    for i in start_positions:
        forward_start[i] = forward_start.get(i, 0) + 1

    json.dump(forward_start, open("data/Rhau/Reverse/start_positions_reverse.json", 'w'), indent=4)

    end_start = {}
    for i in end_positions:
        end_start[i] = end_start.get(i, 0) + 1 
    json.dump(end_start, open("data/Rhau/Reverse/end_positions_reverse.json", 'w'), indent=4)


    reverse_counts = {"R3" : 1407261, "R4" : 1020586, "R5" : 1220302, "R6" : 948704, "R7" : 957320, "R8" : 110, "non_matched" : 230623}
    graph_barcode_counts(reverse_counts, "Reverse Read Barcode Distribution", "data/Rhau/Reverse/barcode_distribution_reverse.png")
    graph_cut_positions(start_positions, end_positions, "data/Rhau/Reverse/cut_positions_distribution_reverse.png")





