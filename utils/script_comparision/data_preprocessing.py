import os
import re
import time

def find_file_pairs(folders):

    file_pairs = []

    current_dir = os.getcwd()

    for folder in folders:
        file_list = os.listdir(os.path.join(current_dir, folder))
        for file in file_list:
            if "R1" in file:
                r1_file = os.path.join(current_dir, folder, file)
                r2_file = r1_file.replace("R1", "R2")
                if os.path.exists(r2_file):
                    file_pairs.append((r1_file, r2_file))

    with open("file_pairs.txt", "w") as f:
        for r1, r2 in file_pairs:
            f.write(f"{r1}\t{r2}\n")

def drop_mismatched_entries(longer_data, shorter_data):
    longer_keys = set(longer_data.keys())
    shorter_keys = set(shorter_data.keys())

    number_to_remove = len(longer_keys) - len(shorter_keys)
    mismatched_keys = list(longer_keys - shorter_keys)

    removed_entries = set()

    for i in range(number_to_remove):
        key_to_remove = mismatched_keys[i]
        removed_entries.add(longer_data[key_to_remove])
        del longer_data[key_to_remove]

    return removed_entries

def check_file_pairs(file_pairs_file, output_dir = "cleaned_pairs"):

    pairs = []
    log = []

    current_dir = os.getcwd()
    output_path = os.path.join(current_dir, output_dir)

    os.makedirs(output_path, exist_ok=True)

    new_pairs = []

    with open(file_pairs_file, "r") as f:
        pairs = [line.strip().split("\t") for line in f]

    for r1, r2 in pairs:

        print(f"Processing pair: {r1} and {r2}")

        start_time = time.time()

        try:
            r1_data, r1_dropped = parseFastqFile(r1)
            r2_data, r2_dropped = parseFastqFile(r2)
        except Exception as e:
            print(f"Error processing files {r1} and {r2}: {e}")
            log.append({
                "r1": r1,
                "r2": r2,
                "r1_entries": 0,
                "r2_entries": 0,
                "r1_dropped": 0,
                "r2_dropped": 0,
                "length_mismatch": 0,
                "time_taken_seconds": 0,
                "error": str(e)})
            continue

        mismatch_length = len(r1_data) - len(r2_data)

        if len(r1_data) != len(r2_data):
            if len(r1_data) > len(r2_data):
                removed_entries = drop_mismatched_entries(r1_data, r2_data)
                r1_dropped.update(removed_entries)
            else:
                removed_entries = drop_mismatched_entries(r2_data, r1_data)
                r2_dropped.update(removed_entries)

        end_time = time.time()

        log.append({
            "r1": r1,
            "r2": r2,
            "r1_entries": len(r1_data) + len(r1_dropped),
            "r2_entries": len(r2_data) + len(r2_dropped),
            "r1_dropped": len(r1_dropped),
            "r2_dropped": len(r2_dropped),
            "length_mismatch": mismatch_length,
            "time_taken_seconds": end_time - start_time
        })

        r1_output_file = os.path.join(output_path, os.path.basename(r1))
        r2_output_file = os.path.join(output_path, os.path.basename(r2))

        with open(r1_output_file, "w") as f:
            for entry in r1_data.values():
                f.write(entry + "\n")
        with open(r2_output_file, "w") as f:
            for entry in r2_data.values():
                f.write(entry + "\n")

        new_pairs.append((r1_output_file, r2_output_file))

    with open("processing_log_2.csv", "w") as f:
        f.write("r1,r2,r1_entries,r2_entries,r1_dropped,r2_dropped,length_mismatch,time_taken_seconds\n")
        for log_entry in log:
            f.write(f"{log_entry['r1']},{log_entry['r2']},{log_entry['r1_entries']},{log_entry['r2_entries']},{log_entry['r1_dropped']},{log_entry['r2_dropped']},{log_entry['length_mismatch']},{log_entry['time_taken_seconds']}\n")

    with open("cleaned_file_pairs_2.txt", "w") as f:
        for r1, r2 in new_pairs:
            f.write(f"{r1}\t{r2}\n")


def read_file_pairs(file_pairs_file):
    pairs = []
    with open(file_pairs_file, "r") as f:
        pairs = [line.strip().split("\t") for line in f]
    return pairs

def parseFastqFile(filePath, flip = False, **kwargs ):

    cull_minlength = kwargs.get("cull_minlength", 8)
    cull_maxlength = kwargs.get("cull_maxlength", 256)
    fileSeperator = r"\+"

    #Load and Read the file
    with open(filePath, "r") as file:
        #Example Entry after regex as a tuple
        #('NB501061:163:HVYLLAFX3:1:11101:1980:1063' - Run ID,  
        # '1:N:0:GTATTATCT+CATATCGTT' - Additional Run Information, 
        # 'TGTAGACTATTCTCACTCTTCTTGTCTGGTTCCTCCGCGTCCGACGTGTGGTGGAGGTTCGGTCGACG', - DNA Sequence 
        # 'AAAAAEE<EE<EEEEEEEEEEAEEEE/EEAEEEEA//AA/EAAEEEEEEEEAEEEEEE/EA<EEEEE6' - DNA Quality Score)
        regexExpression = r"@([A-Z0-9.:\-\s]+)(?:\s)([A-Z0-9:+\/-]+)(?:\s*)([CODONS]{cull_minlength,cull_maxlength})(?:\s+SPLIT\s*)([!-I]{cull_minlength,cull_maxlength})".replace("cull_minlength", str(cull_minlength)).replace("cull_maxlength", str(cull_maxlength))
        regexExpression = regexExpression.replace("CODONS", "ATGCN" )
        
        cleanedFileSeparator = fileSeperator.replace("\\\\", "\\")
        regexExpression = regexExpression.replace("SPLIT", f"[{cleanedFileSeparator}]?")
        entries = re.findall(regexExpression, file.read())

    ent_dict = {}

    dropped_entries = set()

    for entry in entries:

        run_id, additional_info, dna_seq, quality_score = entry

        if len(dna_seq) != len(quality_score):
            dropped_entries.add(entry)
            continue

        ent_dict[run_id] = f"@{run_id} {additional_info}\n{dna_seq}\n+\n{quality_score}"

    return ent_dict, dropped_entries

if __name__ == "__main__":
    folders = ['/home/proxima/Desktop/Side_Projects/P3ANUT/data/Original_Sequencing/M7-Gal', 
               '/home/proxima/Desktop/Side_Projects/P3ANUT/data/Original_Sequencing/M7-Man', 
               "/home/proxima/Desktop/Side_Projects/P3ANUT/data/Original_Sequencing/M7-Nat"]
    find_file_pairs(folders)

    check_file_pairs("missed_files.txt")