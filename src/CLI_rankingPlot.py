"""
CLI and standalone Python function for the Abundance Ranking Plot.

Usage (CLI):
    python CLI_rankingPlot.py \
        --file1 path/to/file1.csv \
        --file2 path/to/file2.csv \
        --output path/to/output.png \
        [--slope 1.0] \
        [--points 100] \
        [--percent-or-count #] \
        [--count-file1] [--count-file2] \
        [--log-scale] \
        [--export-graph]
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")  # non-interactive backend — no GUI required
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

import numpy as np

# Import shared data logic from the GUI module without launching tkinter
sys.path.insert(0, os.path.dirname(__file__))
from rankingPlot import supportingLogic


def run_ranking_plot(
    file1path: str,
    file2path: str,
    output_path: str,
    slope: float = 1.0,
    b: float = 0.0,
    points: float = 100,
    percent_or_count: str = "#",
    count_file1: bool = True,
    count_file2: bool = False,
    export_graph: bool = True,
    log_scale: bool = False,
) -> dict:
    """
    Run the Abundance Ranking Plot logic and optionally save the figure.

    Parameters
    ----------
    file1path       : Path to the first CSV file.
    file2path       : Path to the second CSV file.
    output_path     : Destination path for the exported PNG (required when
                      export_graph=True).
    slope           : Slope of the separation line.
    b               : Y-intercept of the separation line.
    points          : Number/percentage cutoff for sequence selection.
    percent_or_count: "#" for a fixed count or "%" for a percentage.
    count_file1     : Include File 1 sequences in the ranking.
    count_file2     : Include File 2 sequences in the ranking.
    export_graph    : Whether to save the figure to output_path.
    log_scale       : Use a log scale on the Y axis.

    Returns
    -------
    dict with keys:
        "x"           – x-axis ranks
        "y"           – y-axis ranks
        "seqs"        – sequence labels
        "above_count" – sequences above the separation line
        "below_count" – sequences below the separation line
    """
    # ------------------------------------------------------------------ data
    data = supportingLogic.csvComparision(file1path, file2path)
    (x, y), seqs = supportingLogic.gatherScatterData(
        data, count_file1, count_file2, percent_or_count, points
    )
    above_count, below_count = supportingLogic.aboveBelowCounts(x, y, slope, b)

    print(f"Sequences loaded: {len(seqs)}")
    print(f"Above line: {above_count}  |  Below line: {below_count}")

    # --------------------------------------------------------- optional graph
    if export_graph:
        fig, axes = plt.subplots(
            2, 2,
            figsize=(6,4),
            gridspec_kw={"height_ratios": [5, 1], "width_ratios": [1, 5]},
        )
        ax1, ax2 = axes[0, 0], axes[0, 1]
        ax3, ax4 = axes[1, 0], axes[1, 1]

        fig.suptitle("Abundance Ranking Plot", fontsize=14)

        # --- ax1: file name legend ------------------------------------------
        trimmed = [os.path.basename(file1path), os.path.basename(file2path)]
        ax1.set_xlim(0, 3)
        ax1.set_ylim(0, 1)
        ax1.text(1, 0.5, f"File 1 - {trimmed[0]}", ha="center", va="center",
                 fontsize=8, rotation=90)
        ax1.text(2, 0.5, f"File 2 - {trimmed[1]}", ha="center", va="center",
                 fontsize=8, rotation=90)
        ax1.set_yticklabels([])
        ax1.set_yticks([])
        ax1.set_xticklabels([])
        ax1.set_xticks([])
        ax1.set_frame_on(False)
        ax1.set_ylabel("File Names", fontsize=10)

        # --- ax2: main scatter + separation line ----------------------------
        ax2.scatter(x, y, c="black", linewidth=0.5, s=10)

        if len(y) != 0:
            x_line = np.arange(0, min(max(x), round(max(y) / slope)) + 1, 1)
        else:
            x_line = np.arange(0, points, 1)

        y_line = slope * x_line + b

        if len(y) != 0:
            y_end  = min(max(x), round(max(y) / slope)) * slope + b
            max_x  = max(x) if len(x) != 0 else points
            max_y  = max(y) if len(y) != 0 else points
        else:
            y_end  = points * slope + b
            max_x  = points
            max_y  = points

        ax2.plot(x_line, y_line, c="green", linewidth=1)

        stacked = np.stack((x_line, y_line)).T
        patch1 = Polygon(
            [*stacked, [x_line[-1], max_y], [0, max_y]],
            closed=True, color="orange", alpha=0.3,
        )
        patch2 = Polygon(
            [*stacked, [max_x, y_end], [max_x, 0], [0, 0]],
            closed=True, color="lightgray", alpha=0.3,
        )
        ax2.add_patch(patch1)
        ax2.add_patch(patch2)

        ax2.set_xlabel("Sequence Ranking File 1")
        ax2.set_ylabel("Sequence Ranking File 2")
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position("right")

        if log_scale:
            ax2.set_yscale("log")

        # --- ax3: above / below count summary -------------------------------
        corner1 = Polygon([[0, 0], [1, 1], [0, 1]], closed=True,
                           color="orange", alpha=0.5)
        corner2 = Polygon([[0, 0], [1, 1], [1, 0]], closed=True,
                           color="lightgray", alpha=0.5)
        ax3.add_patch(corner1)
        ax3.add_patch(corner2)
        ax3.text(0.25, 0.75, above_count, ha="center", va="center",
                 color="black", fontsize=8)
        ax3.text(0.75, 0.25, below_count, ha="center", va="center",
                 color="black", fontsize=8)
        ax3.set_xlim(0, 1)
        ax3.set_ylim(0, 1)
        ax3.xaxis.set_ticklabels([])
        ax3.get_yaxis().set_visible(False)
        ax3.set_xlabel("Counts")

        # --- ax4: parameter summary -----------------------------------------
        ax4.axis("off")

        fig.tight_layout()
        fig.subplots_adjust(wspace=0, hspace=0)
        fig.savefig(output_path, dpi=100)
        plt.close(fig)
        print(f"Graph saved to: {output_path}")

    return {
        "x": x,
        "y": y,
        "seqs": seqs,
        "above_count": above_count,
        "below_count": below_count,
    }




# --------------------------------------------------------------------------- #
#  CLI entry point                                                             #
# --------------------------------------------------------------------------- #

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="CLI_rankingPlot",
        description="Generate an Abundance Ranking Plot from two CSV files.",
    )
    p.add_argument("--file1",            required=True,
                   help="Path to the first CSV file.")
    p.add_argument("--file2",            required=True,
                   help="Path to the second CSV file.")
    p.add_argument("--output",           required=True,
                   help="Output path for the PNG graph (or CSV results).")
    p.add_argument("--slope",            type=float, default=1.0,
                   help="Slope of the separation line (default: 1.0).")
    p.add_argument("--b",                type=float, default=0.0,
                   help="Y-intercept of the separation line (default: 0.0).")
    p.add_argument("--points",           type=float, default=100,
                   help="Number/percentage cutoff for sequences (default: 100).")
    p.add_argument("--percent-or-count", choices=["%", "#"], default="#",
                   dest="percent_or_count",
                   help="Use a fixed count (#) or percentage (%%) cutoff "
                        "(default: #).")
    p.add_argument("--count-file1",      action="store_true", default=True,
                   help="Include File 1 sequences in the ranking (default: on).")
    p.add_argument("--no-count-file1",   action="store_false", dest="count_file1",
                   help="Exclude File 1 sequences from the ranking.")
    p.add_argument("--count-file2",      action="store_true", default=False,
                   help="Include File 2 sequences in the ranking (default: off).")
    p.add_argument("--log-scale",        action="store_true", default=False,
                   help="Use a log scale on the Y axis.")
    p.add_argument("--export-graph",     action="store_true", default=True,
                   help="Save the figure to --output (default: on).")
    p.add_argument("--no-export-graph",  action="store_false", dest="export_graph",
                   help="Skip saving the figure.")
    return p


def main():
    parser = _build_parser()
    args = parser.parse_args()

    run_ranking_plot(
        file1path       = args.file1,
        file2path       = args.file2,
        output_path     = args.output,
        slope           = args.slope,
        b               = args.b,
        points          = args.points,
        percent_or_count= args.percent_or_count,
        count_file1     = args.count_file1,
        count_file2     = args.count_file2,
        export_graph    = args.export_graph,
        log_scale       = args.log_scale,
    )


if __name__ == "__main__":
    # main()

    slope = 1.0
    b: float = 0.0
    percent_or_count: str = "#"
    count_file1: bool = True
    count_file2: bool = False
    export_graph: bool = True
    log_scale: bool = True

    filePaths = []

    for direction in ["Reverse"]:

        dir_short = "F" if direction == "Forward" else "R"

        for i in [4,5,6,7]:
            for length in ["", "_trimmed"]:
                for count in [10,25,50,100,1000]:

                    input_lambda = lambda run_num: f"data/Rhau/{direction}/R{run_num}_Rhau18_12aa_{dir_short}/R{run_num}_Rhau18_12aa_{dir_short}{length}.csv"
                    output_lambda = lambda run_num: f"data/Rhau/{direction}/R{run_num}_Rhau18_12aa_{dir_short}/C_{count}/"

                    os.makedirs(output_lambda(i), exist_ok=True)

                    file1 = input_lambda(i)
                    file2 = input_lambda(3)
                    output = os.path.join(output_lambda(i), f"rankingPlot_{dir_short}{length}_C{count}.png")

                    #filePaths.append((file1, file2, output, count))


    for i in [4,5,6,7]:
        for length in ["", "_trimmed"]:
            for count in [10,25,50,100,1000]:

                input_lambda = lambda run_num: f"data/Rhau/matched/R{run_num}_Rhau18_12aa/R{run_num}_A30{length}.csv"
                output_lambda = lambda run_num: f"data/Rhau/matched/R{run_num}_Rhau18_12aa/C_{count}/"

                os.makedirs(output_lambda(i), exist_ok=True)

                file1 = input_lambda(i)
                file2 = input_lambda(3)
                output = os.path.join(output_lambda(i), f"rankingPlot_M{length}_C{count}.png")

                filePaths.append((file1, file2, output, count))


    for file1, file2, output, count in filePaths:
        t = run_ranking_plot(
            file1path       = file1,
            file2path       = file2,
            output_path     = output,
            slope           = slope,
            b               = b,
            points          = count,
            percent_or_count= percent_or_count,
            count_file1     = count_file1,
            count_file2     = count_file2,
            export_graph    = export_graph,
            log_scale       = log_scale,
        )
        
        csv_output = output.replace(".png", ".csv")
        with open(csv_output, "w") as f:
            f.write("sequence,rank_file1,rank_file2\n")
            for seq, x_val, y_val in zip(t["seqs"], t["x"], t["y"]):
                f.write(f"{seq},{x_val},{y_val}\n")
        
        pass


