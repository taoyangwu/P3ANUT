"""
CLI and standalone Python function for the Upset Plot.

Usage (CLI):
    python CLI_upsetplot.py \
        --files file1.csv file2.csv file3.csv \
        [--graph-output upset.png] \
        [--export-intersection] \
        [--intersection-output intersection.csv] \
        [--selection true false true]

The --selection flag accepts one boolean (true/false) per file, indicating
which files participate in the intersection query. Defaults to all True.
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from upsetPlot import upsetPlot


def _build_key(selection: list[bool]) -> int:
    """Convert a boolean list to the binary key used by upsetPlot."""
    key = 0
    for i, include in enumerate(reversed(selection)):
        if include:
            key += 2 ** i
    return key


def run_upset_plot(
    files: list[str],
    graph_output: str = None,
    export_intersection: bool = False,
    intersection_output: str = None,
    selection: list[bool] = None,
) -> dict:
    """
    Run the Upset Plot logic programmatically.

    Parameters
    ----------
    files                : List of paths to CSV files (must each have a
                           'sequence' column).
    graph_output         : Optional path to save the figure as a PNG.
    export_intersection  : If True, write the intersection CSV to
                           intersection_output.
    intersection_output  : Path for the intersection CSV output. Required
                           when export_intersection=True.
    selection            : List of booleans (one per file) indicating which
                           files are included in the intersection query.
                           Defaults to all True.

    Returns
    -------
    dict with keys:
        "intersection_counts" – raw intersection count list
        "file_lengths"        – number of sequences in each file
        "total_union_size"    – size of the union across all files
        "key"                 – integer key derived from selection
        "intersection_size"   – number of sequences in the queried intersection
    """
    if len(files) < 2:
        raise ValueError("At least two files are required.")

    if selection is None:
        selection = [True] * len(files)

    if len(selection) != len(files):
        raise ValueError(
            f"selection length ({len(selection)}) must match "
            f"files length ({len(files)})."
        )

    key = _build_key(selection)

    # ------------------------------------------------------------------ data
    intersection_dict = {}
    intersection_counts, file_lengths, total_union_size = upsetPlot.fileSetComparision(
        files, intersection_dict
    )

    intersection_size = intersection_dict.get(key, 0)
    print(f"Files loaded        : {len(files)}")
    print(f"Total union size    : {total_union_size}")
    print(f"Intersection key    : {key}  (selection: {selection})")
    print(f"Intersection size   : {intersection_size}")

    # -------------------------------------------------------- export CSV
    if export_intersection:
        if not intersection_output:
            raise ValueError(
                "intersection_output must be provided when "
                "export_intersection=True."
            )
        upsetPlot.fileOutput(files, intersection_output, key)

    # -------------------------------------------------------- build graph
    if graph_output:
        sorted_axis = np.argsort(
            [t * -1 for t in intersection_counts], kind="stable"
        )
        sorted_counts = np.take(intersection_counts, sorted_axis)

        fig = Figure(figsize=(14, 6), dpi=100)
        fig.suptitle("Upset Plot")
        axes = fig.subplots(
            2, 2,
            gridspec_kw={"height_ratios": [4, 1], "width_ratios": [1, 4]},
        )
        ax1, ax2 = axes[0, 0], axes[0, 1]
        ax3, ax4 = axes[1, 0], axes[1, 1]

        # --- ax2: intersection bar chart -----------------------------------
        max_count = total_union_size
        labels = [
            f"{round(x / max_count * 100)}%" for x in sorted_counts
        ]
        for i, count in enumerate(sorted_counts):
            rect = Rectangle((i, 0), 0.8, count, color="black", alpha=0.5)
            ax2.add_patch(rect)
            ax2.text(
                i + 0.4, count, labels[i],
                ha="center", va="bottom", color="black", fontsize="x-small",
            )
        ax2.set_xlim(0, len(sorted_counts))
        ax2.set_ylim(0, max(sorted_counts) * 1.2 if sorted_counts.any() else 1)
        ax2.set_ylabel("Number of intersections")
        ax2.yaxis.set_label_position("right")
        ax2.yaxis.tick_right()
        ax2.yaxis.label.set_size(10)
        ax2.tick_params(axis="y", labelsize=7)
        ax2.yaxis.set_major_formatter(
            plt.FuncFormatter(lambda x, _: f"{x:.1e}")
        )
        ax2.xaxis.set_visible(False)

        # --- ax4: set association dots + lines -----------------------------
        ax4.set_xticklabels([])
        ax4.set_xticks([])
        ax4.set_xlim(0, len(sorted_axis))

        positives, negatives, connecting_lines = [], [], []
        for i, key_val in enumerate(sorted_axis):
            bits = np.flip(np.power(2, np.arange(len(files))))
            matched = np.sign(np.bitwise_and(key_val, bits))
            connecting_lines.append([])
            for j, v in enumerate(matched):
                if v == 1:
                    positives.append([i + 0.4, j])
                    connecting_lines[-1].append(j)
                else:
                    negatives.append([i + 0.4, j])
        for i, line in enumerate(connecting_lines):
            if len(line) > 1:
                ax4.plot([i + 0.4, i + 0.4], [line[0], line[-1]], c="green")

        file_labels = [os.path.basename(f) for f in files]
        background = np.array([x % 2 == 0 for x in range(len(files))]) * len(sorted_axis)
        ax4.barh(file_labels, background, align="center", color="black", alpha=0.1)
        ax4.get_yaxis().set_visible(False)
        ax4.set_xlabel("Set Association")
        ax4.xaxis.label.set_size(10)

        if positives:
            pos_arr = np.array(positives)
            ax4.scatter(pos_arr[:, 0], pos_arr[:, 1], c="green", s=50)
        if negatives:
            neg_arr = np.array(negatives)
            ax4.scatter(neg_arr[:, 0], neg_arr[:, 1], c="black", alpha=0.5, s=30, zorder=2)
            ax4.scatter(neg_arr[:, 0], neg_arr[:, 1], c="white", alpha=1, s=20, zorder=3)

        # --- ax3: file size horizontal bars --------------------------------
        reference = "ABCDEFGHIJKLMNOP"
        short_names = [reference[i] for i in range(len(files))]
        ax3.barh(short_names, file_lengths, color="black")
        ax3.invert_xaxis()
        ax3.yaxis.tick_right()
        ax3.yaxis.set_label_position("right")
        ax3.set_xlabel("File size")
        ax3.xaxis.label.set_size(10)
        ax3.tick_params(axis="x", labelsize=8)

        # --- ax1: file name legend -----------------------------------------
        trimmed = [os.path.basename(f) for f in files]
        ax1.set_xlim(0, len(files) + 1)
        for i in range(len(files)):
            ax1.text(
                i + 1, 0.5,
                f"{short_names[i]} - {trimmed[i]}",
                ha="center", va="center", fontsize=8, rotation=90,
            )
        ax1.set_yticklabels([])
        ax1.set_yticks([])
        ax1.set_xticklabels([])
        ax1.set_xticks([])
        ax1.set_ylabel("File Names")
        ax1.yaxis.label.set_size(10)

        fig.savefig(graph_output, dpi=100)
        plt.close("all")
        print(f"Graph saved to: {graph_output}")

    return {
        "intersection_counts": intersection_counts,
        "file_lengths": file_lengths,
        "total_union_size": total_union_size,
        "key": key,
        "intersection_size": intersection_size,
    }


# --------------------------------------------------------------------------- #
#  CLI entry point                                                             #
# --------------------------------------------------------------------------- #

def _str_to_bool(value: str) -> bool:
    if value.lower() in ("true", "1", "yes"):
        return True
    if value.lower() in ("false", "0", "no"):
        return False
    raise argparse.ArgumentTypeError(f"Boolean value expected, got: {value}")


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="CLI_upsetplot",
        description="Generate an Upset Plot from multiple CSV files.",
    )
    p.add_argument(
        "--files", required=True, nargs="+",
        help="Paths to the input CSV files (must each contain a 'sequence' column).",
    )
    p.add_argument(
        "--graph-output", default=None, dest="graph_output",
        help="Optional path to save the upset plot figure as a PNG.",
    )
    p.add_argument(
        "--export-intersection", action="store_true", default=False,
        dest="export_intersection",
        help="Export the intersection sequences to a CSV file.",
    )
    p.add_argument(
        "--intersection-output", default=None, dest="intersection_output",
        help="Path for the intersection CSV output "
             "(required when --export-intersection is set).",
    )
    p.add_argument(
        "--selection", nargs="+", type=_str_to_bool, default=None,
        metavar="BOOL",
        help="One true/false per file indicating which files are included "
             "in the intersection query (default: all true). "
             "Example: --selection true false true",
    )
    return p


def main():
    parser = _build_parser()
    args = parser.parse_args()

    run_upset_plot(
        files                = args.files,
        graph_output         = args.graph_output,
        export_intersection  = args.export_intersection,
        intersection_output  = args.intersection_output,
        selection            = args.selection,
    )


if __name__ == "__main__":
    #main()

    #data/Rhau/matched/R7_Rhau18_12aa/C_1000/rankingPlot_M_trimmed_C1000.png
    forward_lambda = lambda I, C, T: f"data/Rhau/Forward/R{I}_Rhau18_12aa_F/C_{C}/rankingPlot_F{T}_C{C}.csv"
    matched_lambda = lambda I, C, T: f"data/Rhau/matched/R{I}_Rhau18_12aa/C_{C}/rankingPlot_M{T}_C{C}.csv"
    reverse_lambda = lambda I, C, T: f"data/Rhau/Reverse/R{I}_Rhau18_12aa_R/C_{C}/rankingPlot_R{T}_C{C}.csv"


    for i in [4,5,6,7]:
        for count in [10, 25, 50, 100,1000]:
            for t in ["", "_trimmed"]:
                file1 = forward_lambda(i, count, t)
                file2 = matched_lambda(i, count, t)
                file3 = reverse_lambda(i, count, t)
                output = f"data/Rhau/UpsetPlots/R{i}_C{count}{t}.png"

                os.makedirs(os.path.dirname(output), exist_ok=True)

                run_upset_plot(
                    files=[file1, file2, file3],
                    graph_output=output,
                    export_intersection=False,
                    intersection_output=None,
                    selection=[True, True, True],
                )

