

"""
Pairwise Levenshtein and Hamming distance calculator.

Given a list of sequences, computes every pair's distances and writes:
  - A CSV with columns: seq_a, seq_b, levenshtein, hamming
  - A PNG confusion-matrix heatmap for each metric
"""

import csv
import os
import itertools

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# --------------------------------------------------------------------------- #
#  Distance functions                                                          #
# --------------------------------------------------------------------------- #

def _levenshtein(a: str, b: str) -> int:
    """Standard dynamic-programming Levenshtein (edit) distance."""
    n, m = len(a), len(b)
    dp = list(range(m + 1))
    for i in range(1, n + 1):
        prev = dp[:]
        dp[0] = i
        for j in range(1, m + 1):
            if a[i - 1] == b[j - 1]:
                dp[j] = prev[j - 1]
            else:
                dp[j] = 1 + min(prev[j - 1], prev[j], dp[j - 1])
    return dp[m]


def _hamming(a: str, b: str) -> int:
    """
    Hamming distance. Sequences of different lengths are padded with a null
    character so the longer sequence counts each extra character as a
    mismatch (distance +1 per extra character).
    """
    max_len = max(len(a), len(b))
    a = a.ljust(max_len, "\x00")
    b = b.ljust(max_len, "\x00")
    return sum(ca != cb for ca, cb in zip(a, b))


# --------------------------------------------------------------------------- #
#  Main analysis function                                                      #
# --------------------------------------------------------------------------- #

def analyse_sequences(
    sequences: list[str],
    csv_output: str = "distances",
    matrix_output_prefix: str = "distance_matrix",
) -> dict:
    """
    Compute pairwise Levenshtein and Hamming similarities for all sequences.

    Similarity is normalized to [0, 1]:
        similarity = 1 - distance / max(len(a), len(b))

    Parameters
    ----------
    sequences             : List of string sequences to compare.
    csv_output            : Path prefix for the similarity matrix CSVs.
                            Two files will be created:
                            <prefix>_levenshtein.csv and <prefix>_hamming.csv
    matrix_output_prefix  : File-path prefix for the heatmap PNGs.
                            Two files will be created:
                            <prefix>_levenshtein.png and <prefix>_hamming.png

    Returns
    -------
    dict with keys:
        "lev_sim_matrix"  – NxN numpy array of Levenshtein similarities
        "ham_sim_matrix"  – NxN numpy array of Hamming similarities
        "lev_matrix"      – NxN numpy array of Levenshtein distances
        "ham_matrix"      – NxN numpy array of Hamming distances
    """
    n = len(sequences)
    if n < 2:
        raise ValueError("At least two sequences are required.")

    lev_matrix     = np.zeros((n, n), dtype=int)
    ham_matrix     = np.zeros((n, n), dtype=int)
    lev_sim_matrix = np.ones((n, n), dtype=float)   # diagonal = 1.0
    ham_sim_matrix = np.ones((n, n), dtype=float)

    for i, j in itertools.combinations(range(n), 2):
        lev = _levenshtein(sequences[i], sequences[j])
        ham = _hamming(sequences[i], sequences[j])
        max_len = max(len(sequences[i]), len(sequences[j]))

        lev_sim = 1.0 - lev / max_len if max_len > 0 else 1.0
        ham_sim = 1.0 - ham / max_len if max_len > 0 else 1.0

        lev_matrix[i, j]     = lev_matrix[j, i]     = lev
        ham_matrix[i, j]     = ham_matrix[j, i]     = ham
        lev_sim_matrix[i, j] = lev_sim_matrix[j, i] = lev_sim
        ham_sim_matrix[i, j] = ham_sim_matrix[j, i] = ham_sim

    # ------------------------------------------------- similarity matrix CSVs
    _write_matrix_csv(lev_sim_matrix, sequences, path=f"{csv_output}_levenshtein.csv")
    _write_matrix_csv(ham_sim_matrix, sequences, path=f"{csv_output}_hamming.csv")

    # --------------------------------------------------------- heatmaps
    _save_heatmap(
        lev_sim_matrix, sequences,
        title="Levenshtein Similarity",
        output_path=f"{matrix_output_prefix}_levenshtein.png",
        value_fmt=".2f",
        colorbar_label="Similarity",
    )
    _save_heatmap(
        ham_sim_matrix, sequences,
        title="Hamming Similarity",
        output_path=f"{matrix_output_prefix}_hamming.png",
        value_fmt=".2f",
        colorbar_label="Similarity",
    )

    return {
        "lev_sim_matrix": lev_sim_matrix,
        "ham_sim_matrix": ham_sim_matrix,
        "lev_matrix": lev_matrix,
        "ham_matrix": ham_matrix,
    }


def _write_matrix_csv(matrix: np.ndarray, sequences: list[str], path: str):
    """Write an NxN similarity matrix CSV with sequences as row/column headers."""
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([""] + list(sequences))        # header row
        for i, seq in enumerate(sequences):
            writer.writerow([seq] + [f"{v:.4f}" for v in matrix[i]])
    print(f"CSV saved to: {path}")


def _save_heatmap(
    matrix: np.ndarray,
    labels: list[str],
    title: str,
    output_path: str,
    value_fmt: str = ".2f",
    colorbar_label: str = "Similarity",
):
    """Render and save a square similarity matrix as a heatmap PNG."""
    n = len(labels)
    fig, ax = plt.subplots(figsize=(max(6, n // 2), max(5, n - 1 // 2)))

    im = ax.imshow(matrix, cmap="viridis", aspect="auto", vmin=0, vmax=1)
    #fig.colorbar(im, ax=ax, label=colorbar_label)

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_title(title)

    # Annotate each cell with the similarity value
    for i in range(n):
        for j in range(n):
            ax.text(
                j, i, format(matrix[i, j], value_fmt),
                ha="center", va="center",
                color="white" if matrix[i, j] < 0.5 else "black",
                fontsize=max(6, 10 - n // 8),
            )

    fig.tight_layout()
    fig.savefig(output_path, dpi=100)
    plt.close(fig)
    print(f"Heatmap saved to: {output_path}")

def _load_sequences(file_path: str, higher_only:bool = True) -> list[str]:
    """Load sequences from a CSV file with a 'seq' column."""
    sequences = []
    with open(file_path, "r") as f:
        reader = csv.reader(f)
        
        for row in reader:

            if higher_only and row[1] <= row[2]:
                sequences.append(row[0])
            elif not higher_only:
                sequences.append(row[0])
    return sequences


# --------------------------------------------------------------------------- #
#  Entry point                                                                 #
# --------------------------------------------------------------------------- #

def main():
    forward_lambda = lambda I, C, T: f"data/Rhau/Forward/R{I}_Rhau18_12aa_F/C_{C}/rankingPlot_F{T}_C{C}.csv"
    matched_lambda = lambda I, C, T: f"data/Rhau/matched/R{I}_Rhau18_12aa/C_{C}/rankingPlot_M{T}_C{C}.csv"
    reverse_lambda = lambda I, C, T: f"data/Rhau/Reverse/R{I}_Rhau18_12aa_R/C_{C}/rankingPlot_R{T}_C{C}.csv"


    for i in [4,5,6,7]:
        for count in [100]:
            for t in ["", "_trimmed"]:
                for dir in ["Forward", "Matched", "Reverse"]:
                    file = forward_lambda(i, count, t) if dir == "Forward" else (matched_lambda(i, count, t) if dir == "Matched" else reverse_lambda(i, count, t))
                
                    output = os.path.join(os.path.dirname(file), f"distance_analysis_{dir}{t}_C{count}")

                    analyse_sequences(
                        sequences=_load_sequences(file, higher_only=True),
                        csv_output=output,
                        matrix_output_prefix=output,
                    )


if __name__ == "__main__":
    main()