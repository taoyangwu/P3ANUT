import argparse
import csv
import os
import re
import shutil
import traceback
import time

import matlab.engine
import json
import re


def parse_forward_files(pairs_file: str) -> list[str]:
    forward_files: list[str] = []
    with open(pairs_file, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            parts = re.split(r"\t+|\s+", line)
            if not parts:
                continue
            forward_files.append((os.path.abspath(os.path.expanduser(parts[0])), False))
    return forward_files

def parse_files(pairs_file: str) -> list[tuple[str, str]]:
    file_pairs: list[tuple[str, str]] = []
    with open(pairs_file, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            parts = re.split(r"\t+|\s+", line)
            if len(parts) < 2:
                continue
            forward_file = os.path.abspath(os.path.expanduser(parts[0]))
            reverse_file = os.path.abspath(os.path.expanduser(parts[1]))
            file_pairs.append((forward_file, False))
            file_pairs.append((reverse_file, True))
    return file_pairs


def pick_largest_text_file(folder: str) -> str:
    candidates = []
    for name in os.listdir(folder):
        full_path = os.path.join(folder, name)
        if os.path.isfile(full_path) and name.lower().endswith(".txt"):
            candidates.append(full_path)

    if not candidates:
        raise FileNotFoundError(f"No text files found in {folder}")

    return max(candidates, key=os.path.getsize)


def count_sequences(good_file: str) -> int:
    total = 0
    with open(good_file, "r") as handle:
        for line in handle:
            tokens = re.split(r"\s+", line.strip())
            if len(tokens) < 2:
                continue
            try:
                total += int(tokens[1])
            except ValueError:
                continue
    return total

def analyze_good_file(good_file: str) -> dict:
    
    tau_score_lambda = lambda seq : "TATTCTCACTCTTCT" in seq and "GGTGGAGGTTCG" in seq
    
    retention_count = 0
    tau_score = 0
    sequence_length_score = 0
    
    upsilion_score = 0
    phi_score = 0
    
    with open(good_file, "r") as handle:
        for line in handle:
            seq = line.strip().split()[0]
            retention_count += 1
            tau_score += 1 if tau_score_lambda(seq) else 0
            sequence_length_score += 1 if len(seq) == 68 else 0
            phi_score += 1 if re.match('[ACGT]*TATTCTCACTCTTCT[ACGT]{27}GGTGGAGGTTCG[ACGT]*', seq) else 0
            upsilion_score += 1 if re.match('[ACGT]{7}TATTCTCACTCTTCT[ACGT]{27}GGTGGAGGTTCG[ACGT]{7}', seq) else 0

    return retention_count, tau_score, sequence_length_score, upsilion_score, phi_score

def _flip_file(file):
    
    flip_bases = {"T": "A", "A": "T", "C": "G", "G": "C", "N": "N"}
    fb_t = str.maketrans(flip_bases)
    filedir, filename = os.path.split(file)
    
    def _load(file_path):

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
            
        return entries
    
    enties = _load(file)
    flipped_entries = []
    for entry in enties:
        run_id, additional_info, sequence, quality = entry
        flipped_sequence = sequence[::-1].translate(fb_t)
        flipped_quality = quality[::-1]
        flipped_entry = f"@{run_id} {additional_info}\n{flipped_sequence}\n+\n{flipped_quality}"
        flipped_entries.append(flipped_entry)
        
    with open(os.path.join(filedir, "flipped_file.fastq"), "w") as file:
        file.write("\n".join(flipped_entries))
        
def _trim_file(file, trim_length_start, trim_length_end):
    
    filedir, filename = os.path.split(file)
    
    def _load(file_path):

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
            
        return entries
    
    enties = _load(file)
    flipped_entries = []
    for entry in enties:
        run_id, additional_info, sequence, quality = entry
        flipped_sequence = sequence[trim_length_start:-trim_length_end ] if trim_length_end > 0 else sequence[trim_length_start:]
        flipped_quality = quality[trim_length_start:-trim_length_end] if trim_length_end > 0 else quality[trim_length_start:]
        flipped_entry = f"@{run_id} {additional_info}\n{flipped_sequence}\n+\n{flipped_quality}"
        flipped_entries.append(flipped_entry)
       
       #os.path.join(filedir, "flipped_file.fastq"), 
    with open(os.path.join(filedir, "trimmed_file.fastq"), "w") as file:
        file.write("\n".join(flipped_entries))
    

def process_file(eng, forward_fastq: str, cleanup: bool, flip = False, trim = 4) -> dict:
    row = {
        "file_fastq": forward_fastq, 
        "sequence_count": "",
        "status": "ok",
        "error": "",
        "time_taken" : "",
        "retention_rate" : "",
        "tau_score" : "",
    }

    filedir, filename = os.path.split(forward_fastq)
    if not os.path.exists(forward_fastq):
        row["status"] = "missing_input"
        row["error"] = "Forward FASTQ does not exist"
        return row
    
    if flip:
        try:
            print(f"Flipping sequences in {filename}")
            _flip_file(os.path.join(filedir, filename))
            filename = "flipped_file.fastq"
        except Exception as exc:
            row["status"] = "error"
            row["error"] = f"Error flipping file: {type(exc).__name__}: {exc}"
            return row
        
    if trim > 0:
        try:
            print(f"Trimming {trim} bases from start and end of sequences in {filename}")
            _trim_file(os.path.join(filedir, filename), trim, 0)
            filename = "trimmed_file.fastq"
        except Exception as exc:
            row["status"] = "error"
            row["error"] = f"Error Trimming file: {type(exc).__name__}: {exc}"
            return row

    bc_folder = os.path.join(filedir, filename.replace(".fastq", "_BC"))

    try:
        shutil.rmtree(bc_folder, ignore_errors=True)
        
        start_time = time.time()

        print(f"Running Step1 on {filename}")
        eng.Step1("inname", filename, "indir", filedir, "indelmut",'on',nargout=0)
        print(f"Finished Step1 for {filename}")

        largest_bc_txt = pick_largest_text_file(bc_folder)
        row["largest_bc_txt"] = largest_bc_txt

        bc_dir, bc_name = os.path.split(largest_bc_txt)
        print(f"Running Step2 on {largest_bc_txt}")
        eng.Step2(
            "inname",
            bc_name,
            "indir",
            bc_dir,
            "start",
            "TCTTGT",
            "end",
            "TTCGAT",
            "uplimit",
            16,
            "DOWNlimit",
            8,
            "fixerr",
            1000,
            "badmax",
            2,
            "ACGTOnly",
            False,
            nargout=0,
        )
        print(f"Finished Step2 for {largest_bc_txt}")
        
        time_taken = time.time() - start_time
        row["time_taken"] = time_taken
        
        row["retention_rate"], row["tau_score"], row["sequence_length_score"], row["upsilion_score"], row["phi_score"] = analyze_good_file(largest_bc_txt)


    except Exception as exc:
        row["status"] = "error"
        row["error"] = f"{type(exc).__name__}: {exc}"

    finally:
        if cleanup:
            shutil.rmtree(bc_folder, ignore_errors=True)

    return row


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run Rebollo Step1/Step2 on forward FASTQ files from filePairs.txt and count resulting sequences."
    )
    parser.add_argument(
        "--pairs",
        default="filePairs.txt",
        help="Path to tab-separated file pairs file (forward in first column).",
    )
    parser.add_argument(
        "--output",
        default="rebollo_forward_counts.csv",
        help="CSV file to write per-forward sequence counts.",
    )
    parser.add_argument(
        "--rebollo-dir",
        default="RebolloScripts",
        help="Path to the Rebollo MATLAB scripts directory.",
    )
    parser.add_argument(
        "--keep-output",
        action="store_true",
        help="Keep generated *_BC output folders instead of deleting them.",
    )
    args = parser.parse_args()

    pairs_file = "filePairs.txt"
    rebollo_dir = "RebolloScripts"
    output_csv = "rebollo_forward_counts_trimmed.csv"

    forward_files = parse_files(pairs_file)
    
    # forward_files = [["data/Full/forclean_PID-1309-GAL-LB1-PC_S93_R1_001.fastq", False],
    #                  ["data/Full/r2ffastq_PID-1309-GAL-LB1-PC_S93_R2_001.fastq", True]]
    
    if not forward_files:
        raise RuntimeError(f"No forward files found in {pairs_file}")

    print(f"Starting MATLAB engine. Files to process: {len(forward_files)}")
    eng = matlab.engine.start_matlab()

    results = []
    try:
        eng.cd(rebollo_dir, nargout=0)

        for idx, (fastq_file, flip) in enumerate(forward_files, start=1):
            
            print(f"[{idx}/{len(forward_files)}] Processing: {fastq_file}")
            row = process_file(eng, fastq_file, cleanup=not args.keep_output, flip=flip, trim=4)
            results.append(row)
            
            with open("rebollo_sanity.json", "w", newline="") as file:
                json.dump(results, file, indent=4)

    except Exception:
        print("Fatal error while running pipeline:")
        traceback.print_exc()
        raise
    finally:
        eng.quit()

    with open("rebollo_sanity.json", "w", newline="") as file:
        json.dump(results, file, indent=4)

    print(f"\nDone. Wrote results to: {output_csv}")


if __name__ == "__main__":
    main()
