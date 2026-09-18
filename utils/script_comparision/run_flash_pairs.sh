#!/usr/bin/env bash

set -u -o pipefail

PAIRS_FILE="${1:-file_pairs.txt}"
MERGED_DIR="${2:-merged}"
FLASH_BIN="${3:-./flash}"

if [[ ! -f "$PAIRS_FILE" ]]; then
  echo "Error: pairs file not found: $PAIRS_FILE" >&2
  exit 1
fi

# Fall back to the common local build path if ./flash is not in the current directory.
if [[ ! -x "$FLASH_BIN" && -x "FLASH-1.2.11/flash" ]]; then
  FLASH_BIN="FLASH-1.2.11/flash"
fi

if [[ ! -x "$FLASH_BIN" ]]; then
  echo "Error: flash executable not found or not executable: $FLASH_BIN" >&2
  echo "Tip: pass the flash path as the 3rd argument, e.g. ./run_flash_pairs.sh file_pairs.txt merged FLASH-1.2.11/flash" >&2
  exit 1
fi

mkdir -p "$MERGED_DIR"

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

  log_file="$MERGED_DIR/${sample}.flash.log"

  echo "[$total] Running FLASH for sample: $sample"
  if "$FLASH_BIN" --max-overlap 75 -d "$MERGED_DIR" -o "$sample" "$r1" "$r2" 2>&1 | tee "$log_file"; then
    ((success++))
  else
    ((failed++))
    echo "FLASH failed for sample: $sample" >&2
  fi
done < "$PAIRS_FILE"

echo
echo "Done. Total processed: $total, Success: $success, Failed: $failed"

if ((failed > 0)); then
  exit 2
fi
