# P3ANUT: Python Pipeline for Peptide Analysis through a Normative Unified Toolset

P3ANUT processes peptide phage display sequencing data, from raw FASTQ reads
through to the comparison plots. The interface is a single page: you drag blocks
onto a canvas, connect them together, and press Run. Every block is also usable
on its own from the command line or from your own Python code.

---

## Contents

- [Installation](#installation)
- [Running the interface](#running-the-interface)
- [Building a pipeline](#building-a-pipeline)
- [Connection types](#connection-types)
- [The blocks](#the-blocks)
- [Saving and reopening a pipeline](#saving-and-reopening-a-pipeline)
- [Configuration](#configuration)
- [Command line usage](#command-line-usage)
- [Repository layout](#repository-layout)

---

## Installation

### Prerequisites

- Conda (Miniconda or Anaconda) for environment management

### 1. Install Conda

Visit the [Miniconda documentation](https://www.anaconda.com/docs/getting-started/miniconda/main)
or [anaconda.com/download](https://www.anaconda.com/download) and install the
version for your operating system.

### 2. Create the environment

On Windows use the "Anaconda Prompt" or "Miniconda Prompt"; on macOS and Linux
any terminal will do.

```bash
conda create -n P3ANUT python=3.12
conda activate P3ANUT
conda install pyyaml numpy pandas matplotlib scikit-learn logomaker anytree conda-forge::python-levenshtein
pip install PyQt6
```

### 3. Verify and activate

```bash
conda env list       # P3ANUT should be listed
conda activate P3ANUT
```

### Alternative: pip

With Python 3.12 available:

```bash
pip install -r requirements.txt
```

---

## Running the interface

From the repository root:

```bash
python P3ANUT.py
```

To open straight into an existing pipeline:

```bash
python P3ANUT.py my_pipeline.p3g.yaml
```

The window has the canvas on the left and, on the right, the Run button, the
block palette, a key to the connection shapes, and the run log.

---

## Building a pipeline

1. **Add blocks.** Drag a block from the palette on the right onto the canvas.
2. **Connect them.** Drag from an output on the right edge of one block to an
   input on the left edge of another. Inputs are always on the left, outputs
   always on the right.
3. **Configure a block.** Click it. A popup lists only the settings that block's
   script actually reads, not the whole configuration file.
4. **Run.** Press Run at the top of the right hand bar.

### Connection rules

- One output can feed **many** inputs; an input accepts **one** connection.
- Two ports can only be joined if their connection types match, so a wrong
  connection cannot be made by accident. If a connection is refused, the reason
  appears in the status bar at the bottom of the window.
- A Ranking Output can only be connected to an Output block.

### Other canvas controls

| Action | How |
| --- | --- |
| Open a block's parameters | Click the block |
| Move a block | Drag it |
| Delete a block or connection | Select it and press Delete |
| Zoom | Ctrl (or Cmd) and scroll |
| Grow or shrink a Run Unifier / Upset Plot | The **+** and **−** buttons on the block |

### Running

Press Run and the progress bar fills as blocks complete, showing
*completed blocks / total blocks* and the name of the block in progress. Output
blocks write their files the moment they execute, so results appear as the run
proceeds rather than all at the end.

**If a block fails, the whole run stops immediately.** Branches that do not
depend on the failed block are not carried on with. The error is shown in a
dialog and the full traceback in the run log, and the run's working directory is
kept so the partial results can be inspected. Pressing Run again starts the
pipeline from the beginning — nothing is resumed or reused from an earlier run.

---

## Connection types

Each connection type has its own shape and colour. The number of sides grows
with the number of pipeline stages the data has passed through, so how far along
a piece of data is can be read at a glance. A two-sided shape is skipped so
nothing is mistaken for a plain line.

| Type | Shape | What travels on it |
| --- | --- | --- |
| File Input | circle | A path to a file you selected |
| Literal Input | circle | An exact value: a number, boolean or short string |
| JSON | triangle | Parsed reads, from the Paired Assembler or the FASTA block |
| Meta Data | triangle | The Paired Assembler's run statistics |
| Sequence Counter | 4 sides | A counted sequence table |
| Unifier | 5 sides | A unified, normalised sequence table |
| Output PNG | 6 sides | A rendered figure |
| Ranking Output | 7 sides | A rank change record, for an Output block only |
| Any | circle | The Output block's ports, which take on whatever is wired in |

---

## The blocks

### File Input
Selects a file from disk.

- **Outputs:** File → *File Input*

### Literal Input
Holds an exact value for blocks that take one as an input, such as Reverse,
Target Length, Top N or Slope. The value's type is chosen in the popup.

- **Outputs:** Value → *Literal Input*

### FASTA
Handles the legacy ingestion path: older sequencing files that have no per-base
quality line, unlike FASTQ's four-lines-per-record format. It reads the file and
converts it into the same JSON representation the rest of the pipeline uses.

Which files per record, whether to reverse complement, and any motif trimming
are all set in the popup; the number of lines per record is what distinguishes a
FASTA-style file (2) from a FASTQ-style one (4).

There is no automatic format detection: use the FASTA block for older,
score-less files and the Paired Assembler for standard FASTQ files.

- **Inputs:** File → *File Input*, Reverse → *Literal Input*
- **Outputs:** Output → *JSON*

### Paired Assembler
Merges paired FASTQ reads into a single corrected set of sequences.

- **Inputs:** Forward File → *File Input*, Reverse File → *File Input*
- **Outputs:** Output Data → *JSON*, Meta Data → *Meta Data*

### Sequence Counter
Produces the final table of counted sequences. A sequence is kept when it

1. is exactly the target number of bases long, set with a Literal Input, and
2. carries both the front and back barcodes, which are set in the popup rather
   than wired in.

Three outputs are produced on **every** run — the counted data, a graph of the
sequence length distribution, and a logo plot. None of them is optional, and
none is an alternative to another.

The length distribution graph shows the lengths actually present in the run's
data, with the target length marked and the share of reads landing on it, and
within one base of it, reported.

The **Surrounding** option folds in sequences that are exactly one base off the
target length, matching each to whichever correctly sized sequence it is closest
to — the most common sequence one base away — instead of dropping it. This uses
the same logic as the FASTA block.

Each Sequence Counter keeps its own barcodes and target length, so two of them
in the same pipeline can be configured differently.

- **Inputs:** Data → *JSON*, Target Length → *Literal Input*
- **Outputs:** Data → *Sequence Counter*, Length Distribution Graph → *Output PNG*,
  Logo Plot → *Output PNG*

### Run Unifier
Merges several counted sequence files into one normalised table. Use the **+**
and **−** buttons on the block to add and remove file inputs; each one is
labelled File 1, File 2 and so on.

- **Inputs:** File 1…N → *Sequence Counter*
- **Outputs:** Output → *Unifier*

### Volcano Plot
Compares two unified runs by ratio and p-value and exports one quadrant.
Quadrant 1 — above both thresholds — is the default.

- **Inputs:** File 1 → *Unifier*, File 2 → *Unifier*, Ratio → *Literal Input*,
  P-Value → *Literal Input*
- **Outputs:** Graph Plot → *Output PNG*, Data → *Unifier*

### Upset Plot
Shows how the sequences of several runs overlap. Takes any number of Run Unifier
outputs, added with the **+** and **−** buttons. By default the exported
intersection is the one across all connected files; the `selection` parameter in
the popup takes one true/false per file to change which are included.

- **Inputs:** File 1…N → *Unifier*
- **Outputs:** Plot → *Output PNG*, Data → *Unifier*

### Ranking Plot
Compares two unified runs by rank, separating them with a line.

- **Inputs:** File A → *Unifier*, File B → *Unifier*, Top N → *Literal Input*,
  Slope (M) → *Literal Input*
- **Outputs:**
  - Above-Line Data → *Unifier* — every sequence ranked above the line. This is
    an ordinary Unifier connection, so it chains onwards into another Volcano
    Plot, Upset Plot or Run Unifier.
  - Graph Plot → *Output PNG*
  - Ranking → *Ranking Output* — a CSV of sequence, File A index and File B
    index, recording the rank change for this step only. It can only be
    connected to an Output block.

The y-intercept, the percent-versus-count cutoff mode, the per-file
include/exclude toggles and the log scale are all in the popup, alongside Top N
and Slope which are wired in.

### Output
Saves whatever it receives to the folder and base file name set on the block,
and passes the same data on downstream, so it can sit in the middle of a chain
as well as at the end.

The file extension is **not** yours to choose — it follows from the type of
connection feeding the block, so a Unifier connection writes `.csv`, an Output
PNG writes `.png`, and so on. Files are written the moment the block executes.

- **Inputs:** Data Input → *any connection*
- **Outputs:** Data Output → *any connection*

---

## Saving and reopening a pipeline

**File → Save pipeline** writes the blocks, their connections and their
configured values to a readable YAML file, conventionally named `*.p3g.yaml`.
**File → Open pipeline** restores it.

This exists so a pipeline does not have to be rebuilt by hand. Reopening a file
gives you the same graph to keep editing and run again; it does not resume a
partly finished run, and it does not reuse any previous run's results.

---

## Configuration

`config.yaml` is the template every block is seeded from. It has a section per
block, and each parameter follows the same pattern:

```yaml
sequenceCount:
  sequenceStart:
    Description: The starting peptide barcode a counted sequence must contain
    Type: String
    UsedFunctions: countSequences
    Value: SHSS
```

`Type` is one of `String`, `Int`, `Float`, `Boolean`, `Choice` (with an
`Options` list) or `BooleanList`.

**The template supplies defaults only.** When a block is placed on the canvas it
takes its own private copy of them. Editing one block's popup never changes
another block of the same type, and never changes `config.yaml`. Two Sequence
Counters in one pipeline can use different barcodes and different target
lengths, and both are stored in the saved pipeline file.

`config_HK.yaml` is the same template with the barcodes, merge tuning and target
length used for the HK datasets.

---

## Command line usage

Every block's underlying function stays runnable from the command line and from
your own code, independently of the interface. Each script also exposes a plain
Python function if you would rather import it.

```bash
# Paired assembly
python src/CLI_pairedAssembler.py -f forward.fastq -r reverse.fastq \
    -o assembled.json -y config.yaml

# Legacy FASTA ingestion (--parse-only does just the step the FASTA block uses)
python utils/FASTA_fileConversion.py --filePath reads.fa --parse-only parsed.json

# Sequence counting: writes the CSV, the length distribution and the logo plot
python src/CLI_sequenceCount.py -f assembled.json -o results/ \
    --target-length 68 -y config.yaml --include-surrounding

# Unify several counted runs
python src/CLI_runUnifier.py -f run1.csv run2.csv -o unified.csv

# Volcano plot
python src/CLI_VolcanoPlot.py -a unifiedA.csv -b unifiedB.csv \
    -o quadrant.csv -g volcano.png -r 1.0 -p 1.3 -q 1

# Upset plot
python src/CLI_upsetplot.py --files a.csv b.csv c.csv \
    --graph-output upset.png --export-intersection \
    --intersection-output intersection.csv --selection true false true

# Ranking plot, including the above-line data and the rank change CSV
python src/CLI_rankingPlot.py --file1 a.csv --file2 b.csv \
    --output ranking.png --points 100 --slope 1.0 \
    --ranking-output ranking.csv --above-output above.csv
```

Add `--help` to any of them for the full set of options.

---

## Repository layout

| Path | Contents |
| --- | --- |
| `P3ANUT.py` | Launcher for the interface |
| `config.yaml` | Configuration template, one section per block |
| `config_HK.yaml` | The template with the HK dataset settings |
| `src/p3anut_ui/` | The PyQt interface: canvas, blocks, graph model and runner |
| `src/` | The pipeline modules and their `CLI_*.py` entry points |
| `utils/` | Shared helpers and the analysis scripts and notebooks |
| `scratch/` | Temporary holding area for data files, images and logs |

`utils/` holds the sequence helpers the pipeline imports — `FASTA_fileConversion.py`
and `visualizationGraphs.py` — alongside the standalone analysis scripts and
notebooks that used to live in `dev_Tools/`.

`scratch/` is a temporary holding location, not a permanent home, for the data
files, images and logs that came out of `dev_Tools/` and the loose files that
were sitting in `utils/`.

---

### Software Development Team

This software application is developed by Ethan Koland, Liam Tucker, Jasmyn
Gooding, and Taoyang Wu.

### Reference

Please cite the associated paper:

Uncovering and Correcting Errors in Peptide Phage Display Library Sequencing

by

Liam Tucker, Ethan Koland, Jasmyn Gooding, Hassan Boudjelal, Derek T. Warren,
David Baker, Maria Marin, Taoyang Wu, Chris J. Morris.
