import argparse
import json
import os
import re
import sys
import time
from pathlib import Path

# Allow importing from project root.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))



BC_FILE_RE = re.compile(r"^BC[1-4]\.txt$", re.IGNORECASE)
RUN_NAME_RE = re.compile(r"^(.*)_R([12])_", re.IGNORECASE)


def _extract_run_name(folder_name):
    match = RUN_NAME_RE.search(folder_name)
    if match:
        return match.group(1), match.group(2)
    return folder_name, None


def _find_largest_bc_file(folder_path):
    largest_path = None
    largest_size = -1

    for entry in os.listdir(folder_path):
        if not BC_FILE_RE.match(entry):
            continue
        candidate = os.path.join(folder_path, entry)
        try:
            size = os.path.getsize(candidate)
        except OSError:
            continue
        if size > largest_size:
            largest_size = size
            largest_path = candidate

    return largest_path


def _iter_bc_sequences(file_path):
    with open(file_path, "r") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            parts = stripped.split()
            if len(parts) < 2:
                continue
            sequence = parts[0].upper()
            quality = parts[-1]
            if len(sequence) != len(quality):
                continue
            yield sequence, quality


def _min_quality_score(quality):
    try:
        return min(ord(char) - 33 for char in quality)
    except ValueError:
        return 0


def _build_patterns(start_barcode, end_barcode):
    upsilion_statement = (
        "[TAGC]{7}start_barcode[TGCA]{27}end_barcode[TAGC]{7}"
        .replace("start_barcode", start_barcode)
        .replace("end_barcode", end_barcode)
    )
    phi_statement = (
        "[TAGC]*start_barcode[TGCA]{27}end_barcode[TAGC]*"
        .replace("start_barcode", start_barcode)
        .replace("end_barcode", end_barcode)
    )
    return re.compile(upsilion_statement), re.compile(phi_statement)


def _reverse_complement(sequence):
    complement = str.maketrans("ACGTN", "TGCAN")
    return sequence.upper().translate(complement)[::-1]


def _barcodes_for_direction(start_barcode, end_barcode, read_direction):
    if read_direction == "2":
        # Reverse reads are reverse-complemented relative to forward reads.
        return _reverse_complement(end_barcode), _reverse_complement(start_barcode)
    return start_barcode, end_barcode


def evaluate_bc_file(
    file_path,
    start_barcode,
    end_barcode,
    target_length,
    upsilion_pattern,
    phi_pattern,
):
    total_reads = 0
    start_counts = 0
    end_counts = 0
    both_counts = 0
    upsilon_score = 0
    phi_score = 0
    target_length_count = 0

    for sequence, quality in _iter_bc_sequences(file_path):
        total_reads += 1
        if "N" in sequence:
            continue
        if _min_quality_score(quality) == 0:
            continue

        if upsilion_pattern.match(sequence):
            upsilon_score += 1
        if phi_pattern.match(sequence):
            phi_score += 1

        contain_start = start_barcode in sequence
        contain_end = end_barcode in sequence

        if contain_start:
            start_counts += 1
        if contain_end:
            end_counts += 1
        if contain_start and contain_end:
            both_counts += 1

        if len(sequence) == target_length:
            target_length_count += 1

    return {
        "total_reads": total_reads,
        "start_counts": start_counts,
        "end_counts": end_counts,
        "both_counts": both_counts,
        "upsilon_score": upsilon_score,
        "phi_score": phi_score,
        "sequence_length_score": target_length_count,
    }


def _normalize_metrics(raw_metrics, file_length, time_taken, source_file):
    denominator = file_length if file_length else 1
    return {
        "time_taken": time_taken,
        "retention_rate": raw_metrics["total_reads"] / denominator,
        "tau_score": raw_metrics["both_counts"] / denominator,
        "upsilon_score": raw_metrics["upsilon_score"] / denominator,
        "phi_score": raw_metrics["phi_score"] / denominator,
        "sequence_length_score": raw_metrics["sequence_length_score"] / denominator,
        "source_bc_file": source_file,
    }


def _average_metrics(forward_metrics, reverse_metrics):
    if not forward_metrics and not reverse_metrics:
        return {}

    keys = [
        "time_taken",
        "retention_rate",
        "tau_score",
        "upsilon_score",
        "phi_score",
        "sequence_length_score",
    ]
    averaged = {}
    for key in keys:
        values = []
        if forward_metrics and forward_metrics.get(key) is not None:
            values.append(forward_metrics[key])
        if reverse_metrics and reverse_metrics.get(key) is not None:
            values.append(reverse_metrics[key])
        averaged[key] = sum(values) / len(values) if values else None

    return averaged


def _load_lengths(length_path):
    with open(length_path, "r") as handle:
        return json.load(handle)


def _load_rebollo_times(rebollo_path):
    if not rebollo_path or not os.path.exists(rebollo_path):
        return {"forward": {}, "reverse": {}}

    with open(rebollo_path, "r") as handle:
        data = json.load(handle)

    if not isinstance(data, list):
        return {"forward": {}, "reverse": {}}

    times = {"forward": {}, "reverse": {}}
    for entry in data:
        if entry.get("status") != "ok":
            continue

        file_fastq = entry.get("file_fastq", "")
        if not file_fastq:
            continue

        try:
            time_taken = float(entry.get("time_taken"))
        except (TypeError, ValueError):
            continue

        run_name, read_direction = _extract_run_name(os.path.basename(file_fastq))
        if not read_direction:
            continue

        if read_direction == "1":
            times["forward"][run_name] = time_taken
        else:
            times["reverse"][run_name] = time_taken

    return times


def build_rebollo_filtered(
    root_folder,
    length_path,
    start_barcode,
    end_barcode,
    target_length,
    rebollo_path=None,
):
    sequence_lengths = _load_lengths(length_path)
    rebollo_times = _load_rebollo_times(rebollo_path)

    forward_data = {}
    reverse_data = {}

    for folder in sorted(os.listdir(root_folder)):
        folder_path = os.path.join(root_folder, folder)
        if not os.path.isdir(folder_path):
            continue

        run_name, read_direction = _extract_run_name(folder)
        if not read_direction:
            continue

        bc_file = _find_largest_bc_file(folder_path)
        if not bc_file:
            continue

        direction_start, direction_end = _barcodes_for_direction(
            start_barcode,
            end_barcode,
            read_direction,
        )
        upsilion_pattern, phi_pattern = _build_patterns(
            direction_start,
            direction_end,
        )

        start_time = time.perf_counter()
        raw_metrics = evaluate_bc_file(
            bc_file,
            direction_start,
            direction_end,
            target_length,
            upsilion_pattern,
            phi_pattern,
        )
        time_taken = time.perf_counter() - start_time

        if read_direction == "1":
            rebollo_time = rebollo_times["forward"].get(run_name)
        else:
            rebollo_time = rebollo_times["reverse"].get(run_name)
        if rebollo_time is not None:
            time_taken += rebollo_time

        length_key = "forward_length" if read_direction == "1" else "reverse_length"
        file_length = sequence_lengths.get(run_name, {}).get(length_key, 0)
        metrics = _normalize_metrics(raw_metrics, file_length, time_taken, bc_file)

        if read_direction == "1":
            forward_data[run_name] = metrics
        else:
            reverse_data[run_name] = metrics

    average_data = {}
    for run_name in sorted(set(forward_data) | set(reverse_data)):
        average_data[run_name] = _average_metrics(
            forward_data.get(run_name),
            reverse_data.get(run_name),
        )

    return {
        "forward": forward_data,
        "reverse": reverse_data,
        "average": average_data,
    }


def main():
    parser = argparse.ArgumentParser(description="Generate rebollo filtered metrics.")
    parser.add_argument(
        "--root-folder",
        default="data/Rebollo_Filtered",
        help="Root folder that contains Rebollo_Filtered subfolders.",
    )
    parser.add_argument(
        "--file-lengths",
        default="file_lengths.json",
        help="Path to file_lengths.json for normalization.",
    )
    parser.add_argument(
        "--rebollo-output",
        default="rebollo_output_ups_phi_trimmed.json",
        help="Path to the original Rebollo output JSON for runtime aggregation.",
    )
    parser.add_argument(
        "--output",
        default="rebollo_filtered.json",
        help="Output JSON file path.",
    )
    parser.add_argument(
        "--start-barcode",
        default="TATTCTCACTCTTCT",
        help="Start barcode sequence.",
    )
    parser.add_argument(
        "--end-barcode",
        default="GGTGGAGGTTCG",
        help="End barcode sequence.",
    )
    parser.add_argument(
        "--target-length",
        type=int,
        default=68,
        help="Target sequence length for sequence_length_score.",
    )

    args = parser.parse_args()

    data = build_rebollo_filtered(
        args.root_folder,
        args.file_lengths,
        args.start_barcode,
        args.end_barcode,
        args.target_length,
        args.rebollo_output,
    )

    output_path = Path(args.output)
    with output_path.open("w") as handle:
        json.dump(data, handle, indent=4)


if __name__ == "__main__":
    main()