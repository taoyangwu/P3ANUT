"""
CLI and standalone Python function for the Sequence Counter.

Usage (CLI):
    python CLI_sequenceCount.py \
        --file pairedAssembleOutput.json \
        --output results/ \
        --target-length 48 \
        [--yamlFile config.yaml] \
        [--include-surrounding] \
        [--sequence-start SHSS] [--sequence-end GGGS] \
        [--method direct] [--encoding ONEHOT]

The block produces exactly three outputs on every run: the counted sequence
CSV, a graph of the sequence length distribution present in the data, and a
logo plot.
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")  # non-interactive backend - no GUI required

import yaml

sys.path.insert(0, os.path.dirname(__file__))
import sequenceCounter


def run_sequence_count(
    data_path: str,
    target_length: int,
    output_directory: str,
    base_name: str = None,
    **config,
) -> dict:
    """
    Run the Sequence Counter logic programmatically.

    Parameters
    ----------
    data_path        : Path to the JSON produced by the Paired Assembler or the
                       FASTA block, or a CSV/list of sequences.
    target_length    : Number of bases a read must contain to be counted.
    output_directory : Folder the three outputs are written to.
    base_name        : Stem for the output file names. Defaults to the input
                       file name.
    **config         : Any sequenceCount configuration value, for example
                       sequenceStart, sequenceEnd, includeSurrounding,
                       minimumCount, method or encoding.

    Returns
    -------
    dict with keys "counts", "dataPath", "lengthPlotPath", "logoPlotPath"
    and "totalReads".
    """
    os.makedirs(output_directory, exist_ok=True)

    if base_name is None:
        base_name = os.path.splitext(os.path.basename(data_path))[0]

    return sequenceCounter.countSequences(
        data_path,
        target_length,
        outputCSV=os.path.join(output_directory, f"{base_name}.csv"),
        lengthPlotPath=os.path.join(output_directory, f"{base_name}_lengthDistribution.png"),
        logoPlotPath=os.path.join(output_directory, f"{base_name}_logoPlot.png"),
        **config,
    )


# --------------------------------------------------------------------------- #
#  Argument validation                                                         #
# --------------------------------------------------------------------------- #

def fileInput(string):
    """Validate that the input exists and is a format the counter can read."""
    if not os.path.isfile(string):
        raise FileNotFoundError(string)

    if string.split(".")[-1].lower() not in ["json", "csv", "fastq", "fq"]:
        raise argparse.ArgumentTypeError("File must be a json, csv or fastq file")

    return string


def dir_path(string):
    """Validate that the path is a directory, creating it when it is missing."""
    os.makedirs(string, exist_ok=True)
    return string


def yaml_path(string):
    """Validate that the path is a loadable YAML file."""
    if not os.path.isfile(string):
        raise FileNotFoundError(string)

    with open(string, 'r') as stream:
        try:
            yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
            raise

    return string


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog='CLI_sequenceCount',
        description='Count the sequences of a target length that carry both barcodes.',
    )
    p.add_argument('-f', '--file', dest='file', metavar='file', type=fileInput,
                   required=True, help='The JSON or CSV file to count')
    p.add_argument('-o', '--output', dest='output', metavar='output', type=dir_path,
                   required=True, help='Folder the three outputs are written to')
    p.add_argument('-t', '--target-length', dest='target_length', type=int,
                   default=None,
                   help='Number of bases a read must contain. Falls back to the '
                        'targetLength value in the YAML file.')
    p.add_argument('--base-name', dest='base_name', default=None,
                   help='Stem for the output file names (default: input file name)')

    p.add_argument('-y', '--yamlFile', metavar='yamlFile', type=yaml_path,
                   help='YAML file holding the sequenceCount defaults')

    #Overrides for the most commonly changed settings
    p.add_argument('--sequence-start', dest='sequenceStart', default=None,
                   help='The starting peptide barcode')
    p.add_argument('--sequence-end', dest='sequenceEnd', default=None,
                   help='The ending peptide barcode')
    p.add_argument('--include-surrounding', dest='includeSurrounding',
                   action='store_true', default=None,
                   help='Fold in sequences one base off the target length')
    p.add_argument('--no-barcodes', dest='matchBarcodes', action='store_false',
                   default=None, help='Count every sequence of the target length')
    p.add_argument('--minimum-count', dest='minimumCount', type=int, default=None,
                   help='Drop sequences seen fewer than this many times')
    p.add_argument('--normalize', dest='normalizeCount', action='store_true',
                   default=None, help='Divide every count by the total read count')
    p.add_argument('--encoding', metavar='encoding', choices=sequenceCounter.validEncodings(),
                   default=None, help='Encoding used by the clustering methods')
    p.add_argument('--method', metavar='method', choices=sequenceCounter.validMethods(),
                   default=None, help='Counting metric to apply')

    p.add_argument('--version', action='version', version='%(prog)s 2.0')
    return p


def main():
    args = _build_parser().parse_args()

    config = {}
    if args.yamlFile:
        with open(args.yamlFile, 'r') as stream:
            yamlSettings = yaml.safe_load(stream)

        config = {k: v.get("Value")
                  for k, v in yamlSettings.get("sequenceCount", {}).items()
                  if v.get("Value") is not None}

    #Anything given on the command line wins over the YAML defaults
    for key in ("sequenceStart", "sequenceEnd", "includeSurrounding", "matchBarcodes",
                "minimumCount", "normalizeCount", "encoding", "method"):
        value = getattr(args, key)
        if value is not None:
            config[key] = value

    target_length = args.target_length or config.get("targetLength")
    if target_length is None:
        raise SystemExit("A target length is required: pass --target-length or "
                         "supply a YAML file with sequenceCount.targetLength set.")

    config.pop("targetLength", None)

    run_sequence_count(
        data_path        = args.file,
        target_length    = target_length,
        output_directory = args.output,
        base_name        = args.base_name,
        **config,
    )


if __name__ == "__main__":
    main()
