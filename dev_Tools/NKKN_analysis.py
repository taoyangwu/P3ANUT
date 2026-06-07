import re

INCLUDE_N = True  # Set to True to include 'N' in the DNA sequences, False to exclude
CULL_MINLENGTH = 10  # Minimum length of DNA sequences to be extracted
CULL_MAXLENGTH = 1000  # Maximum length of DNA sequences to be extracted
FILE_SEPARATOR = "+"  # Separator used in the input file, can be space, tab, etc.
K = {"G","T"}

def extract_sequences(file_path):
    with open(file_path, "r") as file:
        #Example Entry after regex as a tuple
        #('NB501061:163:HVYLLAFX3:1:11101:1980:1063' - Run ID,  
        # '1:N:0:GTATTATCT+CATATCGTT' - Additional Run Information, 
        # 'TGTAGACTATTCTCACTCTTCTTGTCTGGTTCCTCCGCGTCCGACGTGTGGTGGAGGTTCGGTCGACG', - DNA Sequence 
        # 'AAAAAEE<EE<EEEEEEEEEEAEEEE/EEAEEEEA//AA/EAAEEEEEEEEAEEEEEE/EA<EEEEE6' - DNA Quality Score)
        regexExpression = r"@([A-Z0-9.:\-\s]+)(?:\s)([A-Z0-9:+\/-]+)(?:\s*)([CODONS]{cull_minlength,cull_maxlength})(?:\s+SPLIT\s*)([!-I]{cull_minlength,cull_maxlength})".replace("cull_minlength", str(CULL_MINLENGTH)).replace("cull_maxlength", str(CULL_MAXLENGTH))
        regexExpression = regexExpression.replace("CODONS", "ATGCN" if INCLUDE_N else "ATGC")
        
        cleanedFileSeparator = FILE_SEPARATOR.replace("\\\\", "\\")
        regexExpression = regexExpression.replace("SPLIT", f"[{cleanedFileSeparator}]?")
        entries = re.findall(regexExpression, file.read())
        
        
    return [x[2] for x in entries]  # Return only the DNA sequences


def main():
    file_path = "data/MON_BSA/PID-1309-M7-MON-BSA-1_S87_R1_001.fastq"  # Replace with your actual file path
    sequences = extract_sequences(file_path)
    
    nkkn_counts = {}
    dropped_N_count = 0
    dropped_K_count = 0
    total_count = 0
    
    count_of_bases_in_position = [{"A": 0, "T": 0, "G": 0, "C": 0, "N": 0} for _ in range(4)]
    
    
    for seq in sequences:
        first_four = seq[:4]
        
        total_count += 1
        
        if("N" in first_four):
            dropped_N_count += 1
            continue
        
        #Check K
        if first_four[1] not in K or first_four[2] not in K:
            dropped_K_count += 1
            continue
        
        nkkn_counts[first_four] = nkkn_counts.get(first_four, 0) + 1
        
        for i, base in enumerate(first_four):
            if base in count_of_bases_in_position[i]:
                count_of_bases_in_position[i][base] += 1
        
    for nkkn, count in nkkn_counts.items():
        print(f"{nkkn}: {count}")
        
    max_nkkn = max(nkkn_counts, key=nkkn_counts.get)
    min_nkkn = min(nkkn_counts, key=nkkn_counts.get)
    
    print(f"Most common NKKN: {max_nkkn} with count {nkkn_counts[max_nkkn]}")
    print(f"Least common NKKN: {min_nkkn} with count {nkkn_counts[min_nkkn]}")
        
    #Output the count of bases in each position as a csv
    print("Position, A, T, G, C, N")
    for i, counts in enumerate(count_of_bases_in_position):
        print(f"{i+1}, {counts['A']}, {counts['T']}, {counts['G']}, {counts['C']}, {counts['N']}")
        
    pass
        
        
        
if __name__ == "__main__":
    main()