import re
import numpy as np
import random
from matplotlib import pyplot as plt

INCLUDE_N = True  # Set to True to include 'N' in the DNA sequences, False to exclude
CULL_MINLENGTH = 10  # Minimum length of DNA sequences to be extracted
CULL_MAXLENGTH = 1000  # Maximum length of DNA sequences to be extracted
FILE_SEPARATOR = "+"  # Separator used in the input file, can be space, tab, etc.
K = {"G","T"}
K_REV = {"C","A"}
N = {"A","T","G","C"}

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

def _validate_umi(umi):

        validate_lambda = lambda u, T_K: u[1] in T_K and u[2] in T_K and all(base in "ATGC" for base in u)

        if len(umi) == 8:
            return validate_lambda(umi[:4], K) and validate_lambda(umi[4:], K_REV)
        if len(umi) == 4:
            return validate_lambda(umi, K)
        
        raise ValueError(f"UMI length must be either 4 or 8, got {len(umi)}")
    
def _validate_region(region, target_length=68):
        return len(region) == target_length - 47  and all(base in "ATGC" for base in region)

def umi_analysis(sequences, only_front_umi : bool =True, target_length : int = 68):

    umi_fixed = {}
    region_fixed = {}

    umi_lambda = lambda seq: seq[:4] if only_front_umi else seq[:4] + seq[-4:]
    region_lambda = lambda seq: seq[25:-22]
    
    for seq in sequences:
        umi = umi_lambda(seq)
        region = region_lambda(seq)

        t = _validate_umi(umi) 
        t2 = _validate_region(region)

        if not _validate_umi(umi) or not _validate_region(region):
            continue

        umi_fixed.setdefault(umi, {})
        umi_fixed[umi][region] = umi_fixed[umi].get(region, 0) + 1

        region_fixed.setdefault(region, {})
        region_fixed[region][umi] = region_fixed[region].get(umi, 0) + 1

    '''
    #Calculate the mean and std of the number of UMIs per region
    umi_unique_region_counts = [len(regions) for regions in umi_fixed.values()]

    mean_umis_per_region = sum(umi_unique_region_counts) / len(umi_unique_region_counts)
    std_umis_per_region = (sum((x - mean_umis_per_region) ** 2 for x in umi_unique_region_counts) / len(umi_unique_region_counts)) ** 0.5

    print(f"Mean UMIs per region: {mean_umis_per_region}")
    print(f"Standard Deviation of UMIs per region: {std_umis_per_region}")

    region_unique_umi_counts = [len(umis) for umis in region_fixed.values()]
    mean_regions_per_umi = sum(region_unique_umi_counts) / len(region_unique_umi_counts)
    std_regions_per_umi = (sum((x - mean_regions_per_umi) ** 2 for x in region_unique_umi_counts) / len(region_unique_umi_counts)) ** 0.5
    print(f"Mean regions per UMI: {mean_regions_per_umi}")
    print(f"Standard Deviation of regions per UMI: {std_regions_per_umi}")
    '''

    t = set(sequences)
        
    umi_unique_region_counts = np.array([len(regions) for regions in umi_fixed.values()])
    mean_umis_per_region = np.mean(umi_unique_region_counts)
    std_umis_per_region = np.std(umi_unique_region_counts)

    region_unique_umi_counts = np.array([len(umis) for umis in region_fixed.values()])
    mean_regions_per_umi = np.mean(region_unique_umi_counts)
    std_regions_per_umi = np.std(region_unique_umi_counts)

    print(f"Mean UMIs per region: {mean_umis_per_region}")
    print(f"Standard Deviation of UMIs per region: {std_umis_per_region}")
    print(f"Mean regions per UMI: {mean_regions_per_umi}")
    print(f"Standard Deviation of regions per UMI: {std_regions_per_umi}")


    return umi_fixed, region_fixed

def _umi_gt_analysis(sequences, only_front_umi : bool =True, target_length : int = 68):

    def _generate_random_umi(only_front_umi):

        temp_N = list(N)
        temp_K = list(K)
        temp_K_REV = list(K_REV)

        if only_front_umi:
            b1 = random.choice(temp_N)
            b2 = random.choice(temp_K)
            b3 = random.choice(temp_K)
            b4 = random.choice(temp_N)
            return b1 + b2 + b3 + b4
        else:
            b5 = random.choice(temp_N)
            b6 = random.choice(temp_K_REV)
            b7 = random.choice(temp_K_REV)
            b8 = random.choice(temp_N)
            return _generate_random_umi(True) + b5 + b6 + b7 + b8

    region_counts = {}

    for seq in sequences:
        umi = seq[:4] if only_front_umi else seq[:4] + seq[-4:]
        region = seq[25:-22]

        if not _validate_region(region, target_length):
            continue

        region_counts[region] = region_counts.get(region, 0) + 1

    synth_data = [(_generate_random_umi(only_front_umi), seq) for seq in region_counts.keys()]

    umi_fixed = {}
    region_fixed = {}

    for umi, region in synth_data:
        umi_fixed.setdefault(umi, {})
        umi_fixed[umi][region] = umi_fixed[umi].get(region, 0) + 1

        region_fixed.setdefault(region, {})
        region_fixed[region][umi] = region_fixed[region].get(umi, 0) + 1

    umi_unique_region_counts = np.array([len(regions) for regions in umi_fixed.values()])
    mean_umis_per_region = np.mean(umi_unique_region_counts)
    std_umis_per_region = np.std(umi_unique_region_counts)

    region_unique_umi_counts = np.array([len(umis) for umis in region_fixed.values()])
    mean_regions_per_umi = np.mean(region_unique_umi_counts)
    std_regions_per_umi = np.std(region_unique_umi_counts)

    print(f"Mean UMIs per region: {mean_umis_per_region}")
    print(f"Standard Deviation of UMIs per region: {std_umis_per_region}")
    print(f"Mean regions per UMI: {mean_regions_per_umi}")
    print(f"Standard Deviation of regions per UMI: {std_regions_per_umi}")


    return umi_fixed, region_fixed
    pass

def plot_region_distribution(region_fixed):
    counts_per_top_umi = []

    for region, umis in region_fixed.items():
        sorted_umis = sorted(umis.items(), key=lambda x: x[1], reverse=True)
        total_umis = sum(umis.values())

        if total_umis <= 100:
            continue

        for i, (umi, count) in enumerate(sorted_umis):

            if i > 100:
                break

            if len(counts_per_top_umi) <= i:
                counts_per_top_umi.append([])
            counts_per_top_umi[i].append(count/total_umis)

    plt.figure(figsize=(10, 6))
    plt.boxplot(counts_per_top_umi)
    plt.ylabel('Count of Regions Associated with UMI')
    plt.title('Distribution of Region Counts for Top UMIs')
    plt.show()

def NKKN_analysis(sequences):
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



def main():
    file_path = "data/Original_Sequencing/M7-Gal/PID-1309-GAL-BSA-1-PC_S107_R1_001.fastq"  # Replace with your actual file path
    sequences = extract_sequences(file_path)

    umi_fixed, region_fixed = umi_analysis(sequences, only_front_umi=False, target_length=68)

    plot_region_distribution(region_fixed)

    #umi_counts = _umi_gt_analysis(sequences, only_front_umi=False, target_length=68)
    
    
        
        
        
if __name__ == "__main__":
    main()