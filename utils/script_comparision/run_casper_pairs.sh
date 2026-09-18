#!/usr/bin/env bash

set -u -o pipefail

PAIRS_FILE="${1:-file_pairs.txt}"
OUTPUT_DIR="${2:-casper_merged}"
CASPER_BIN="${3:-./casper}"

if [[ ! -f "$PAIRS_FILE" ]]; then
  echo "Error: pairs file not found: $PAIRS_FILE" >&2
  exit 1
fi

if [[ ! -x "$CASPER_BIN" ]]; then
  echo "Error: CASPER executable not found or not executable: $CASPER_BIN" >&2
  echo "Tip: pass the CASPER path as the 3rd argument, e.g. ./run_casper_pairs.sh cleaned_file_pairs.txt casper_merged ./casper" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Resolve paths so the command still works when running CASPER inside OUTPUT_DIR.
CASPER_BIN_ABS="$(readlink -f "$CASPER_BIN")"
OUTPUT_DIR_ABS="$(readlink -f "$OUTPUT_DIR")"

total=0
success=0
failed=0

while IFS= read -r line || [[ -n "$line" ]]; do
  # Skip comments and blank lines.
  [[ -z "${line//[[:space:]]/}" ]] && continue
  [[ "${line#\#}" != "$line" ]] && continue

  r1=""
  r2=""
  read -r r1 r2 _ <<< "$line"

  if [[ -z "$r1" || -z "$r2" ]]; then
    echo "Skipping malformed line: $line" >&2
    ((failed++))
    continue
  fi

  if [[ ! -f "$r1" || ! -f "$r2" ]]; then
    echo "Skipping missing files:" >&2
    echo "  R1: $r1" >&2
    echo "  R2: $r2" >&2
    ((failed++))
    continue
  fi

  ((total++))

  sample="$(basename "$r1")"
  sample="${sample%_R1_001.fastq}"
  sample="${sample%_R1.fastq}"
  sample="${sample%_R1_001.fq}"
  sample="${sample%_R1.fq}"

  log_file="$OUTPUT_DIR_ABS/${sample}.casper.log"

  echo "[$total] Running CASPER for sample: $sample"
  if (
    cd "$OUTPUT_DIR_ABS" &&
    "$CASPER_BIN_ABS" "$r1" "$r2" -o "$sample"
  ) 2>&1 | tee "$log_file"; then
    ((success++))
  else
    ((failed++))
    echo "CASPER failed for sample: $sample" >&2
  fi
done < "$PAIRS_FILE"

echo
echo "Done. Total processed: $total, Success: $success, Failed: $failed"

if ((failed > 0)); then
  exit 2
fi
