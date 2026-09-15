"""
CLI and standalone Python function for the Volcano Plot.

Usage (CLI):
    python CLI_VolcanoPlot.py \
        --InputA fileA.csv \
        --InputB fileB.csv \
        --output volcanoPlot.csv \
        [--graph-output volcano.png] \
        [--ratio 1.0] [--pvalue 1.3] [--quadrant 1]

Quadrant 1, the default, is everything above both the ratio and the p-value
threshold.
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")  # non-interactive backend - no GUI required
import matplotlib.pyplot as plt

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from volcanoPlot import supportingLogic


#Shading applied to each quadrant, indexed by quadrant number
_QUADRANT_COLORS = {1: "tab:red", 2: "tab:blue", 3: "lightgray", 4: "tab:orange"}


def run_volcano_plot(
    file_a: str,
    file_b: str,
    output_path: str = None,
    graph_output: str = None,
    ratio: float = 1.0,
    pvalue: float = 1.3,
    quadrant: int = 1,
    figure_width: float = 8.0,
    figure_height: float = 6.0,
    dpi: int = 300,
) -> dict:
    """
    Run the Volcano Plot logic and optionally save the figure and the
    selected quadrant.

    Parameters
    ----------
    file_a        : Path to the first unified CSV.
    file_b        : Path to the second unified CSV.
    output_path   : Where to write the CSV of the selected quadrant.
    graph_output  : Where to save the figure as a PNG.
    ratio         : A/B ratio threshold, the vertical split.
    pvalue        : -log10(p-value) threshold, the horizontal split.
    quadrant      : Which quadrant is exported on the data output (default 1).
    figure_width  : Figure width in inches.
    figure_height : Figure height in inches.
    dpi           : Resolution of the exported figure.

    Returns
    -------
    dict with keys:
        "data"            - the joined ratio / p-value DataFrame
        "quadrant_data"   - the rows falling in the selected quadrant
        "quadrant_counts" - counts for quadrants 1 to 4
        "output_path"     - the CSV written, when one was requested
        "graph_output"    - the PNG written, when one was requested
    """
    data = supportingLogic.csvComparision(file_a, file_b)

    quadrant_data = supportingLogic.returnQuadrant(data, ratio, pvalue, quadrant)
    quadrant_counts = supportingLogic.quarantCounts(data, ratio, pvalue)

    print(f"Sequences compared : {len(data)}")
    print(f"Quadrant counts    : Q1={quadrant_counts[0]} Q2={quadrant_counts[1]} "
          f"Q3={quadrant_counts[2]} Q4={quadrant_counts[3]}")
    print(f"Exporting quadrant : {quadrant} ({len(quadrant_data)} sequences)")

    if output_path:
        export = supportingLogic.exportDF(file_a, file_b, quadrant_data)
        export.to_csv(output_path)
        print(f"Quadrant data saved to: {output_path}")

    if graph_output:
        x = data['AvB_Ratio'].values
        y = data['-log10(P-Value)'].values

        fig, ax = plt.subplots(figsize=(figure_width, figure_height))

        ax.scatter(x, y, c="black", s=6, alpha=0.4, linewidth=0)

        #Highlight whichever quadrant is being exported
        if len(quadrant_data):
            ax.scatter(quadrant_data['AvB_Ratio'].values,
                       quadrant_data['-log10(P-Value)'].values,
                       c=_QUADRANT_COLORS.get(quadrant, "tab:red"), s=10,
                       linewidth=0, label=f"Quadrant {quadrant}")

        ax.axvline(ratio, color="green", linestyle="--", linewidth=1,
                   label=f"Ratio = {ratio}")
        ax.axhline(pvalue, color="blue", linestyle="--", linewidth=1,
                   label=f"-log10(p) = {pvalue}")

        #Label each quadrant with how many sequences landed in it
        x_max = np.nanmax(x) if len(x) else 1
        y_max = np.nanmax(y) if len(y) else 1
        corners = {
            1: (ratio + (x_max - ratio) / 2, pvalue + (y_max - pvalue) / 2),
            2: (ratio / 2, pvalue + (y_max - pvalue) / 2),
            3: (ratio / 2, pvalue / 2),
            4: (ratio + (x_max - ratio) / 2, pvalue / 2),
        }
        for index, (cx, cy) in corners.items():
            ax.text(cx, cy, f"Q{index}\n{quadrant_counts[index - 1]}",
                    ha="center", va="center", fontsize=9, alpha=0.7)

        ax.set_xlabel("File A / File B ratio")
        ax.set_ylabel("-log10(P-Value)")
        ax.set_title("Volcano Plot")
        ax.legend(fontsize=8)

        fig.tight_layout()
        fig.savefig(graph_output, dpi=dpi)
        plt.close(fig)
        print(f"Graph saved to: {graph_output}")

    return {
        "data": data,
        "quadrant_data": quadrant_data,
        "quadrant_counts": quadrant_counts,
        "output_path": output_path,
        "graph_output": graph_output,
    }


# --------------------------------------------------------------------------- #
#  CLI entry point                                                             #
# --------------------------------------------------------------------------- #

def fileInput(string):
    """Validate that the input is an existing CSV file."""
    if not os.path.isfile(string):
        raise FileNotFoundError(string)

    if string.split(".")[-1].lower() != "csv":
        raise argparse.ArgumentTypeError("File must be a csv file")

    return string


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog='CLI_VolcanoPlot',
        description='Compare two unified CSV files and export a volcano plot quadrant.',
    )
    p.add_argument('-a', "--InputA", metavar='fileA', type=fileInput, required=True,
                   help='The A file to be used in the comparison')
    p.add_argument('-b', "--InputB", metavar='fileB', type=fileInput, required=True,
                   help='The B file to be used in the comparison')
    p.add_argument('-r', "--ratio", metavar='ratio', type=float, default=1.0,
                   help='The A/B ratio threshold (default: 1.0)')
    p.add_argument('-p', "--pvalue", metavar='pvalue', type=float, default=1.3,
                   help='The -log10(p-value) threshold (default: 1.3)')
    p.add_argument('-q', "--quadrant", metavar='quadrant', type=int, default=1,
                   choices=[1, 2, 3, 4], help='Quadrant to export (default: 1)')
    p.add_argument('-o', "--output", metavar='output', type=str,
                   default="volcanoPlot.csv", help='Path for the quadrant CSV')
    p.add_argument('-g', "--graph-output", dest='graph_output', default=None,
                   help='Optional path to save the figure as a PNG')
    p.add_argument('--dpi', type=int, default=300,
                   help='Resolution of the exported figure (default: 300)')
    return p


def main():
    args = _build_parser().parse_args()

    run_volcano_plot(
        file_a       = args.InputA,
        file_b       = args.InputB,
        output_path  = args.output,
        graph_output = args.graph_output,
        ratio        = args.ratio,
        pvalue       = args.pvalue,
        quadrant     = args.quadrant,
        dpi          = args.dpi,
    )


if __name__ == "__main__":
    main()
