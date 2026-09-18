#!/usr/bin/env python3

import argparse
import json
import re
from pathlib import Path


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


    


def main() -> int:
    

    flash_dir = "merged"
    caspar_dir = "casper_merged"
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

            run_name = log_path.name[: -len(".flash.log")]
            summary.setdefault(run_name, {})["flash"] = parse_flash_log(log_path)

            if measure_tau:

                output_file = log_path.replace_suffix(".flash.log").with_suffix("extendedFrags.fastq")

                tau_score = calculate_tau_score(output_file)
                run_name = log_path.name[: -len(".flash.log")]
                summary.setdefault(run_name, {})["tau_score"] = tau_score


            
            flash_count += 1
    else:
        print(f"Warning: FLASH log directory not found, skipping: {flash_dir}")

    if casper_dir.is_dir():
        for log_path in sorted(casper_dir.glob("*.casper.log")):
            run_name = log_path.name[: -len(".casper.log")]
            summary.setdefault(run_name, {})["casper"] = parse_casper_log(log_path)
            casper_count += 1
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
    raise SystemExit(main())
