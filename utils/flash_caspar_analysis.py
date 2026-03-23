#!/usr/bin/env python3

import argparse
import json
import re
from pathlib import Path
import os

from evaluation_comparision import evalutate_fastq_file


TOTAL_RE = re.compile(r"Total pairs:\s*([0-9,]+)")
COMBINED_RE = re.compile(r"Combined pairs:\s*([0-9,]+)")
UNCOMBINED_RE = re.compile(r"Uncombined pairs:\s*([0-9,]+)")
PERCENT_RE = re.compile(r"Percent combined:\s*([0-9]+(?:\.[0-9]+)?)%")
SECONDS_RE = re.compile(r"([0-9]+(?:\.[0-9]+)?)\s+seconds elapsed")

CASPER_TOTAL_RE = re.compile(r"Total number of reads\s*:\s*([0-9,]+)")
CASPER_MERGED_RE = re.compile(
    r"Number of merged reads\s*:\s*([0-9,]+)\s*\(([0-9]+(?:\.[0-9]+)?)%\)"
)
CASPER_UNMERGED_RE = re.compile(
    r"Number of unmerged reads\s*:\s*([0-9,]+)\s*\(([0-9]+(?:\.[0-9]+)?)%\)"
)
CASPER_TIME_RE = re.compile(
    r"TIME for total processing\s*:\s*([0-9]+(?:\.[0-9]+)?)\s*sec"
)


def _to_int(value: str) -> int:
    return int(value.replace(",", ""))


def parse_flash_log(log_path: Path) -> dict:
    metrics = {
        "total_pairs": None,
        "combined_pairs": None,
        "uncombined_pairs": None,
        "percent_combined": None,
        "time_seconds": None,
    }

    for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if metrics["total_pairs"] is None:
            m = TOTAL_RE.search(line)
            if m:
                metrics["total_pairs"] = _to_int(m.group(1))
                continue

        if metrics["combined_pairs"] is None:
            m = COMBINED_RE.search(line)
            if m:
                metrics["combined_pairs"] = _to_int(m.group(1))
                continue

        if metrics["uncombined_pairs"] is None:
            m = UNCOMBINED_RE.search(line)
            if m:
                metrics["uncombined_pairs"] = _to_int(m.group(1))
                continue

        if metrics["percent_combined"] is None:
            m = PERCENT_RE.search(line)
            if m:
                metrics["percent_combined"] = float(m.group(1))
                continue

        m = SECONDS_RE.search(line)
        if m:
            metrics["time_seconds"] = float(m.group(1))

    return metrics


def parse_casper_log(log_path: Path) -> dict:
    metrics = {
        "total_pairs": None,
        "combined_pairs": None,
        "uncombined_pairs": None,
        "percent_combined": None,
        "time_seconds": None,
    }

    for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if metrics["total_pairs"] is None:
            m = CASPER_TOTAL_RE.search(line)
            if m:
                metrics["total_pairs"] = _to_int(m.group(1))
                continue

        if metrics["combined_pairs"] is None:
            m = CASPER_MERGED_RE.search(line)
            if m:
                metrics["combined_pairs"] = _to_int(m.group(1))
                metrics["percent_combined"] = float(m.group(2))
                continue

        if metrics["uncombined_pairs"] is None:
            m = CASPER_UNMERGED_RE.search(line)
            if m:
                metrics["uncombined_pairs"] = _to_int(m.group(1))
                continue

        m = CASPER_TIME_RE.search(line)
        if m:
            metrics["time_seconds"] = float(m.group(1))

    return metrics

def calculate_tau_score(file):
    return evalutate_fastq_file(file, start_barcode="TATTCTCACTCTTCT", end_barcode="GGTGGAGGTTCG")

def calculate_sequence_length(file_path, target_length = 68):

    cull_minlength = 8
    cull_maxlength = 256
    fileSeperator = r"\+"

    #Load and Read the file
    with open(file_path, "r") as file:
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

    sequence_lengths = [0] * (cull_maxlength + 1)

    for entry in entries:
        sequence = entry[2]
        sequence_lengths[len(sequence)] += 1

    return sequence_lengths[target_length] / sum(sequence_lengths)


def main() -> int:
    

    flash_dir = "/home/proxima/Desktop/Side_Projects/FLASH_CASPAR/merged"
    caspar_dir = "/home/proxima/Desktop/Side_Projects/FLASH_CASPAR/casper_merged"
    output_path = "metrics.json"

    flash_dir = Path(flash_dir)
    casper_dir = Path(caspar_dir)
    output_path = Path(output_path)
    measure_tau = True

    summary = {}
    flash_count = 0
    casper_count = 0

    if flash_dir.is_dir():
        for log_path in sorted(flash_dir.glob("*.flash.log")):

            print(f"Parsing FLASH log: {log_path}")

            run_name = log_path.name[: -len(".flash.log")]
            summary.setdefault(run_name, {})["flash"] = parse_flash_log(log_path)

            if measure_tau:

                output_file = os.path.join(log_path.parent, log_path.stem.replace(".flash", "") + ".extendedFrags.fastq")

                

                try:
                    tau_score = calculate_tau_score(output_file)
                    summary[run_name]["flash"] = summary[run_name]["flash"] | tau_score
                    summary[run_name]["flash"]["sequence_length_score"] = calculate_sequence_length(output_file)
                except Exception as e:
                    print(f"Error calculating tau score for {output_file}: {e}")
                    summary[run_name]["flash"]["tau_score"] = None
            
            flash_count += 1
    else:
        print(f"Warning: FLASH log directory not found, skipping: {flash_dir}")

    if casper_dir.is_dir():
        for log_path in sorted(casper_dir.glob("*.casper.log")):

            print(f"Parsing CASPER log: {log_path}")

            run_name = log_path.name[: -len(".casper.log")]
            summary.setdefault(run_name, {})["casper"] = parse_casper_log(log_path)

            

            if measure_tau:

                output_file = os.path.join(log_path.parent, log_path.stem.replace(".casper", "") + ".fastq")

                

                try:
                    tau_score = calculate_tau_score(output_file)
                    summary[run_name]["casper"] = summary[run_name]["casper"] | tau_score
                    summary[run_name]["casper"]["sequence_length_score"] = calculate_sequence_length(output_file)
                except Exception as e:
                    print(f"Error calculating tau score for {output_file}: {e}")
                    summary[run_name]["casper"]["tau_score"] = None

            
            flash_count += 1
    else:
        print(f"Warning: CASPER log directory not found, skipping: {casper_dir}")

    if flash_count == 0 and casper_count == 0:
        raise SystemExit("No .flash.log or .casper.log files were found.")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")

    print(f"Parsed FLASH logs: {flash_count}")
    print(f"Parsed CASPER logs: {casper_count}")
    print(f"Total unique runs in output: {len(summary)}")
    print(f"Wrote metrics to {output_path}")
    return 0


if __name__ == "__main__":
    # from evaluation_comparision import t

    # t()

    raise SystemExit(main())
