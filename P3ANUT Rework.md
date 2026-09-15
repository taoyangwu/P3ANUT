Refactor the code of the implementation of the following multiple stage UI for processing FASTQ (and legacy FASTA) files.

The new UI should use the PyQt framework rather than the current Tkinter framework, replacing `unifiedGUI.py`, `GUI_pairedAssembler.py`, `GUI_sequenceCount.py`, `configPopupV2.py`, `runUnifier.py`, `volcanoPlot.py`, `upsetPlot.py`, and `rankingPlot.py`'s Tkinter UI code, and should exist as one single UI page.

The new UI should utilize a codeless UI to connect the components of each step to the other. A bar on the right should allow the user to drag in blocks to the UI. The left contains the canvas where the user connects the components.

On the top of the right section there should be a Run button. Clicking it runs every block in the graph while displaying a progress bar. Progress is calculated as (number of blocks completed) / (total number of blocks in the graph). If any block errors, stop the entire run immediately — do not continue processing unaffected branches. Partial/resumable execution of the pipeline is out of scope for this version.

Different connection types should be represented by different shapes and colors to prevent the user from making the wrong connections. The number of sides of a connection's shape increases with the number of pipeline stages the data has passed through (a raw File Input connection is a circle; the output of the Paired Assembler and FASTA blocks is a triangle; a 2-sided shape is skipped to avoid confusion with a plain line). Colors should follow the matplotlib "tab"/Tableau color cycle.

One output should be able to connect to multiple inputs, but an input should only accept one connection.

Inputs should sit on the left of each block; outputs should sit on the right.

When a block is clicked, a popup should show the user the parameters that can be configured for that block. These are the same parameters located within the config YAML, but rather than showing every parameter in the config, only the ones actually used within that block's underlying script/functions should be shown. Each block instance stores its own copy of these parameter values, seeded from the YAML template's defaults — editing one block's popup does not affect any other instance of the same block type elsewhere in the graph.

Allow the user to save and load the current block graph (blocks, their connections, and their configured values) in a user-readable file format. The purpose is purely to avoid having to rebuild the pipeline from scratch — reopening a saved file restores the full graph so the user can keep editing and then re-run it; it does not imply resuming a partially-completed run or reusing any previous run's cached output.

## Blocks

### File Input Block
An input block that allows the user to select a file.

### Output Block
A block that saves the data it receives to the file type, folder, and base filename specified within the block, and also acts as a pass-through that forwards the data downstream. The block writes its file to disk as soon as it executes within the run (not batched until the end of the run).

The block should have a folder icon that lets the user select the destination folder and base file name. The file extension is not user-defined — it is set automatically based on the connection type feeding the block.

#### Inputs
- Data Input -> any connection

#### Outputs
- Data Output -> any connection

### Literal Input Block
A block that allows the user to enter an exact value (used for things like booleans, numbers, and short strings that other blocks otherwise take in as inputs — e.g. Reverse, Target Length, Top N, Slope).

### FASTA Block
Handles the legacy FASTA ingestion path: older sequencing files that do not include a per-base quality/read-score line, unlike FASTQ's 4-line-per-record format. This is only the step where the FASTA file is loaded and transformed into the same JSON-style data representation the rest of the pipeline uses.

This block partially calls `utils/FASTA_fileConversion.py` (confirmed as the correct script). That file has many components in it (parsing, amino conversion, logo plotting, comparison scatter plots); this block should only invoke the file-reading/line-parsing part of that script — notably the handling of the different possible read-line layouts (e.g. how many lines make up one entry) — not the amino conversion, counting, or plotting logic.

Routing between this block and the Paired Assembler block is manual: the user is responsible for dragging in the correct block for their file (FASTA Block for older, score-less files; Paired Assembler for standard FASTQ files). There is no automatic file-format detection.

#### Inputs
- File -> File Input connection
- Reverse -> Literal Input connection

#### Outputs
- Output -> JSON connection

### Paired Assembler
Represents the paired assembler function (`multiprocessedPairAssembler.py`). The block should look at its inputs and parse them according to the paired assembler script.

#### Inputs
- Forward File -> File Input connection
- Reverse File -> File Input connection

#### Outputs
- Output Data -> JSON connection
- Meta Data -> Output connection

### Sequence Counter
Rewritten to output only the final CSV of sequences that (a) match the target number of bases, set via a Literal Input block, and (b) match the front and back barcodes, which are configured within the block's parameter popup rather than wired in (i.e. `sequenceStart`/`sequenceEnd` stay config-popup values, not connections). These popup values are scoped per block instance, not shared globally — two Sequence Counter blocks in the same graph can be configured with different target lengths and barcodes. Drop the current multi-file outputs such as `unmatched.csv`, `purged.csv`, and the per-length CSVs.

The logo plot is a third, standard output — Data, Length Distribution Graph, and Logo Plot are all produced on every run, not a togglable alternative.

Add a further output that graphs the distribution of sequence lengths seen in the data, using the `graph_sequence_length_distribution` figure logic in `utils/visualizationGraphs.py`. This function needs to be reworked, not just re-parameterized: today it hardcodes "target length 68" and indexes a single fixed bin — it needs to instead plot the actual sequence-length distribution present in a given run's data, so the user can see how lengths are distributed rather than a single hardcoded value.

A "Surrounding" sequence is a sequence whose length is exactly one base off the target length. The "Surrounding" config option folds these into the output by matching them to whichever correctly-sized sequence they are closest to (the most common sequence one base away), rather than dropping them, using the same logic as the FASTA block.

#### Inputs
- Data -> JSON connection
- Target Length (number of bases) -> Literal Input connection

#### Outputs
- Data -> Sequence Counter connection
- Length Distribution Graph -> Output PNG connection
- Logo Plot -> Output PNG connection

### Run Unifier
This block should have an additional feature: clicking it adds another file input, with a matching subtract button to remove excess connections. The title for each input should be "File N" where N is that input's position.

#### Inputs
- File 1...N -> Sequence Counter connection

#### Outputs
- Output -> Unifier connection

### Volcano Plot
The default behavior should default to quadrant 1.

#### Inputs
- File 1 -> Unifier connection
- File 2 -> Unifier connection
- Ratio -> Literal Input connection
- P-Value -> Literal Input connection

#### Outputs
- Graph Plot -> Output PNG connection
- Data -> Unifier connection

### Upset Plot
Accepts any number of Run Unifier outputs as file inputs. Defaults to the intersection across all connected files; which files are included in the exported intersection should be configurable (mirroring the existing per-file true/false selection already supported by `CLI_upsetplot.py`) via a new `upsetPlot` section added to the config YAML (no such section exists today).

#### Inputs
- File 1...N -> Unifier connection

#### Outputs
- Plot -> Output PNG connection
- Data -> Unifier connection (the data for the selected intersection)

### Ranking Plot
New block comparing two Run Unifier outputs by rank.

#### Inputs
- File A -> Unifier connection
- File B -> Unifier connection
- Top N -> Literal Input connection
- Slope (M) -> Literal Input connection

#### Outputs
- Above-Line Data -> Unifier connection (all sequences ranked above the line — chains into further blocks like any other Unifier connection, e.g. into another Volcano Plot, Upset Plot, or Run Unifier)
- Graph Plot -> Output PNG connection
- Ranking -> Ranking Output connection (a CSV with columns: sequence, File A index, File B index). This output only reports the rank change for this step. The UI should only allow a Ranking Output connection to be wired into an Output block — it cannot be connected into any other block type.

## Additional Features
- Maintain the CLI and programmatic interface to each step, so every block's underlying function stays runnable from code and from the command line, independent of the UI.
- Update `README.md` to document the new PyQt UI, the new/changed blocks (including the new Ranking Plot block and the rewritten Sequence Counter output), and the reorganized `utils/`/`scratch/` folders.
- Replace `configV32.yaml` (and regenerate `config_HK.yaml` to match) with a new YAML config template containing sections for every block, including new sections for Run Unifier, Volcano Plot, Upset Plot, and Ranking Plot (none of these have a config section today — only `global`, `pairedAssembler`, and `sequenceCount` exist in `configV32.yaml`), following the same `Description`/`Type`/`UsedFunctions`/`Value` pattern already used for the existing sections. This template supplies each parameter's default value only — each block instance in a saved graph keeps its own independent copy of the actual values (see above), so different instances of the same block type can diverge from the template defaults without affecting each other. Ranking Plot's full existing CLI surface (slope, y-intercept, the percent-vs-count cutoff mode, the per-file include/exclude toggles, log-scale) should all be exposed as config-popup parameters via this new section, alongside Top N and Slope as the block's wired Literal Inputs.
- Move all Python scripts and Jupyter notebooks (`.py`/`.ipynb`) currently in `dev_Tools/` into `utils/`, updating any hardcoded paths that referenced the old `dev_Tools/` location. Where a script of the same name already exists in both folders (`visualizationGraphs.py`), the `utils/` version is the one to keep.
- Move all other non-Python, non-notebook files (data files, images, logs, etc.) currently in `dev_Tools/` — and any already sitting loose in `utils/`, such as `log(1).txt` and `results.csv` — into a new `scratch/` folder. `scratch/` is only a temporary holding location for this material, not a permanent home.
