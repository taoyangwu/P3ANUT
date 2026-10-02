"""
Block definitions for the P3ANUT node graph.

Every block declares the ports it exposes, the configuration parameters its
underlying script actually reads, and how to execute it. The execute functions
call straight into the pipeline modules, so a block and its command line
equivalent always run the same code.

Data is handed between blocks as a Payload: a value plus the connection type it
travelled on. For everything except a literal, that value is a path to a file in
the run's working directory, which keeps the underlying scripts - most of which
take file paths - usable without modification.
"""

import json
import os
import shutil
import sys
from collections import namedtuple
from dataclasses import dataclass, field
from typing import Callable

from . import connections as C

#The pipeline modules live in src/, one level up from this package, and the
#shared sequence helpers live in utils/ at the repository root.
_SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_REPO_ROOT = os.path.dirname(_SRC)

for _location in (_SRC, _REPO_ROOT):
    if _location not in sys.path:
        sys.path.insert(0, _location)


Payload = namedtuple("Payload", "value type")


@dataclass(frozen=True)
class PortSpec:
    key: str
    label: str
    type: C.ConnectionType


@dataclass
class BlockDefinition:
    key: str
    label: str
    configSection: str
    #The parameters this block's script actually reads. None means the whole
    #section, which is the right answer when the section exists for this block
    #alone.
    configKeys: list = None
    inputs: Callable = field(default=lambda node: [])
    outputs: Callable = field(default=lambda node: [])
    execute: Callable = None
    #Blocks whose input count the user grows and shrinks on the canvas
    dynamicInputs: bool = False
    dynamicLabel: str = "File"
    description: str = ""


def value(inputs, key, default=None):
    """Read a wired input's value, falling back when nothing is connected."""
    payload = inputs.get(key)
    return default if payload is None else payload.value


def _require(inputs, key, blockLabel, portLabel):
    """Read an input that the block cannot run without."""
    payload = inputs.get(key)
    if payload is None or payload.value in (None, ""):
        raise ValueError(f"{blockLabel}: the '{portLabel}' input is not connected")

    return payload.value


# --------------------------------------------------------------------------- #
#  File Input                                                                  #
# --------------------------------------------------------------------------- #

def _executeFileInput(node, inputs, ctx):
    path = node.params.get("filePath")

    if not path:
        raise ValueError("File Input: no file has been selected")

    if not os.path.isfile(path):
        raise FileNotFoundError(f"File Input: {path}")

    return {"file": path}


FILE_INPUT_BLOCK = BlockDefinition(
    key="fileInput",
    label="File Input",
    configSection="fileInput",
    outputs=lambda node: [PortSpec("file", "File", C.FILE_INPUT)],
    execute=_executeFileInput,
    description="Selects a file from disk and feeds it into the pipeline.",
)


# --------------------------------------------------------------------------- #
#  Output                                                                      #
# --------------------------------------------------------------------------- #

def _uniquePath(path, overwrite):
    """Return path, or path with a numeric suffix when it is already taken."""
    if overwrite or not os.path.exists(path):
        return path

    stem, extension = os.path.splitext(path)
    index = 1
    while os.path.exists(f"{stem}_{index}{extension}"):
        index += 1

    return f"{stem}_{index}{extension}"


def _executeOutput(node, inputs, ctx):
    payload = inputs.get("data")

    if payload is None:
        raise ValueError("Output: nothing is connected to the data input")

    directory = node.params.get("outputDirectory") or ""
    baseName = node.params.get("baseFileName") or "output"

    if not directory:
        raise ValueError("Output: no destination folder has been chosen. Use the "
                         "folder button on the block to pick one.")

    os.makedirs(directory, exist_ok=True)

    #The extension is not the user's to choose - it follows from the type of
    #connection feeding the block.
    extension = C.extensionFor(payload.type, payload.value)
    destination = _uniquePath(os.path.join(directory, baseName + extension),
                              node.params.get("overwrite", True))

    if payload.type is C.LITERAL:
        with open(destination, "w") as handle:
            handle.write(str(payload.value))
    else:
        shutil.copyfile(payload.value, destination)

    ctx.log(f"Output: wrote {destination}")

    #Also a pass-through, so the same data carries on downstream
    return {"data": payload.value}


OUTPUT_BLOCK = BlockDefinition(
    key="output",
    label="Output",
    configSection="outputBlock",
    inputs=lambda node: [PortSpec("data", "Data Input", C.ANY)],
    outputs=lambda node: [PortSpec("data", "Data Output", C.ANY)],
    execute=_executeOutput,
    description="Saves whatever it receives to the chosen folder and passes it on. "
                "The file extension follows the connection type.",
)


# --------------------------------------------------------------------------- #
#  Literal Input                                                               #
# --------------------------------------------------------------------------- #

def _executeLiteral(node, inputs, ctx):
    from .config import coerce

    declared = node.params.get("valueType", "String")
    raw = node.params.get("value", "")

    try:
        return {"value": coerce(raw, declared)}
    except ValueError as exc:
        raise ValueError(f"Literal Input: '{raw}' is not a valid {declared}") from exc


LITERAL_BLOCK = BlockDefinition(
    key="literal",
    label="Literal Input",
    configSection="literalInput",
    outputs=lambda node: [PortSpec("value", "Value", C.LITERAL)],
    execute=_executeLiteral,
    description="Holds an exact value - a number, a boolean or a short string - "
                "for blocks that take one as an input.",
)


# --------------------------------------------------------------------------- #
#  FASTA                                                                       #
# --------------------------------------------------------------------------- #

def _executeFasta(node, inputs, ctx):
    from utils.FASTA_fileConversion import fastaToPipelineData

    filePath = _require(inputs, "file", "FASTA", "File")

    params = dict(node.params)
    #A wired Reverse literal overrides the value held in the popup
    reverse = value(inputs, "reverse")
    if reverse is not None:
        params["flip_sequence"] = bool(reverse)

    data = fastaToPipelineData(
        filePath,
        dnaDatatag=params.get("dnaDatatag", "sequences"),
        dataEnteryLength=params.get("dataEnteryLength", 2),
        flip_sequence=params.get("flip_sequence", False),
        motrif_presequence_length=params.get("motrif_presequence_length", 0),
        motrif_postsequence_length=params.get("motrif_postsequence_length", 0),
    )

    destination = ctx.path(node, "fasta.json")
    with open(destination, "w") as handle:
        json.dump(data, handle)

    ctx.log(f"FASTA: parsed {len(data)} reads from {os.path.basename(filePath)}")

    return {"output": destination}


FASTA_BLOCK = BlockDefinition(
    key="fasta",
    label="FASTA",
    configSection="fastaBlock",
    inputs=lambda node: [
        PortSpec("file", "File", C.FILE_INPUT),
        PortSpec("reverse", "Reverse", C.LITERAL),
    ],
    outputs=lambda node: [PortSpec("output", "Output", C.JSON)],
    execute=_executeFasta,
    description="Reads a legacy, score-less sequencing file and converts it into the "
                "JSON representation the rest of the pipeline uses.",
)


# --------------------------------------------------------------------------- #
#  Paired Assembler                                                            #
# --------------------------------------------------------------------------- #

def _executePairedAssembler(node, inputs, ctx):
    import multiprocessedPairAssembler as pairedassembler

    forward = _require(inputs, "forward", "Paired Assembler", "Forward File")
    reverse = value(inputs, "reverse")

    params = dict(node.params)

    #These two are addressed through globalPatameters rather than kwargs, and
    #are kept inside the run directory so a run never litters the cwd.
    globalParameters = {
        "logFile": ctx.path(node, params.pop("logFileName", "logfile.txt")),
        "errorFile": ctx.path(node, params.pop("errorFileName", "errorData.json")),
    }

    data = {}
    metadata = pairedassembler.parse(forward, reverse, data,
                                     globalPatameters=globalParameters,
                                     addToLogFile=True, **params)

    dataPath = ctx.path(node, "pairedAssembly.json")
    with open(dataPath, "w") as handle:
        json.dump(data, handle)

    metaPath = ctx.path(node, "pairedAssemblyMeta.json")
    with open(metaPath, "w") as handle:
        json.dump(metadata, handle, indent=4)

    ctx.log(f"Paired Assembler: {metadata.get('finalCount', len(data))} merged reads")

    return {"data": dataPath, "meta": metaPath}


PAIRED_ASSEMBLER_BLOCK = BlockDefinition(
    key="pairedAssembler",
    label="Paired Assembler",
    configSection="pairedAssembler",
    inputs=lambda node: [
        PortSpec("forward", "Forward File", C.FILE_INPUT),
        PortSpec("reverse", "Reverse File", C.FILE_INPUT),
    ],
    outputs=lambda node: [
        PortSpec("data", "Output Data", C.JSON),
        PortSpec("meta", "Meta Data", C.METADATA),
    ],
    execute=_executePairedAssembler,
    description="Merges paired FASTQ reads into a single corrected set of sequences.",
)


# --------------------------------------------------------------------------- #
#  Sequence Counter                                                            #
# --------------------------------------------------------------------------- #

def _executeSequenceCounter(node, inputs, ctx):
    import sequenceCounter

    dataPath = _require(inputs, "data", "Sequence Counter", "Data")

    params = dict(node.params)
    #A wired Target Length literal overrides the value held in the popup
    targetLength = value(inputs, "targetLength", params.get("targetLength"))

    if targetLength is None:
        raise ValueError("Sequence Counter: no target length. Connect a Literal Input "
                         "to Target Length or set it in the block's parameters.")

    params.pop("targetLength", None)

    #Named after the block, so that a Run Unifier downstream ends up with a
    #distinguishable column per source rather than several called the same thing.
    outputCSV = ctx.path(node, f"{node.id}.csv")
    lengthPlot = ctx.path(node, f"{node.id}_lengthDistribution.png")
    logoPlot = ctx.path(node, f"{node.id}_logoPlot.png")

    sequenceCounter.countSequences(dataPath, targetLength,
                                   outputCSV=outputCSV,
                                   lengthPlotPath=lengthPlot,
                                   logoPlotPath=logoPlot,
                                   **params)

    return {"data": outputCSV, "lengthGraph": lengthPlot, "logoPlot": logoPlot}


SEQUENCE_COUNTER_BLOCK = BlockDefinition(
    key="sequenceCounter",
    label="Sequence Counter",
    configSection="sequenceCount",
    inputs=lambda node: [
        PortSpec("data", "Data", C.JSON),
        PortSpec("targetLength", "Target Length", C.LITERAL),
    ],
    outputs=lambda node: [
        PortSpec("data", "Data", C.SEQUENCE_COUNT),
        PortSpec("lengthGraph", "Length Distribution Graph", C.OUTPUT_PNG),
        PortSpec("logoPlot", "Logo Plot", C.OUTPUT_PNG),
    ],
    execute=_executeSequenceCounter,
    description="Counts the sequences of the target length that carry both barcodes, "
                "and graphs the length distribution and the logo plot.",
)


# --------------------------------------------------------------------------- #
#  Run Unifier                                                                 #
# --------------------------------------------------------------------------- #

def _dynamicInputs(node, connectionType, label):
    return [PortSpec(f"file{i + 1}", f"{label} {i + 1}", connectionType)
            for i in range(max(1, node.inputCount))]


def _collectDynamic(node, inputs, blockLabel):
    """The wired dynamic inputs, in port order, rejecting any gaps."""
    paths = []
    for index in range(max(1, node.inputCount)):
        key = f"file{index + 1}"
        payload = inputs.get(key)
        if payload is None:
            raise ValueError(f"{blockLabel}: input {index + 1} is not connected")
        paths.append(payload.value)

    return paths


def _executeRunUnifier(node, inputs, ctx):
    from runUnifier import merge

    paths = _collectDynamic(node, inputs, "Run Unifier")

    #Named after the block so that a downstream Upset Plot gets a
    #distinguishable column per source rather than several called the same thing.
    destination = ctx.path(node, f"{node.id}.csv")
    merge(paths, destination)

    ctx.log(f"Run Unifier: merged {len(paths)} files")

    return {"output": destination}


RUN_UNIFIER_BLOCK = BlockDefinition(
    key="runUnifier",
    label="Run Unifier",
    configSection="runUnifier",
    #inputCount seeds the starting port count but is driven by the + and -
    #buttons on the block, so it is not offered in the popup as well.
    configKeys=[],
    inputs=lambda node: _dynamicInputs(node, C.SEQUENCE_COUNT, "File"),
    outputs=lambda node: [PortSpec("output", "Output", C.UNIFIER)],
    execute=_executeRunUnifier,
    dynamicInputs=True,
    description="Merges several counted sequence files into one normalised table.",
)


# --------------------------------------------------------------------------- #
#  Volcano Plot                                                                #
# --------------------------------------------------------------------------- #

def _executeVolcano(node, inputs, ctx):
    from CLI_VolcanoPlot import run_volcano_plot

    fileA = _require(inputs, "file1", "Volcano Plot", "File 1")
    fileB = _require(inputs, "file2", "Volcano Plot", "File 2")

    params = dict(node.params)
    ratio = value(inputs, "ratio", params.get("ratio", 1.0))
    pvalue = value(inputs, "pvalue", params.get("pvalue", 1.3))

    graphPath = ctx.path(node, "volcano.png")
    dataPath = ctx.path(node, "volcanoQuadrant.csv")

    run_volcano_plot(
        fileA, fileB,
        output_path=dataPath,
        graph_output=graphPath,
        ratio=float(ratio),
        pvalue=float(pvalue),
        quadrant=int(params.get("quadrant", 1)),
        figure_width=float(params.get("figureWidth", 8.0)),
        figure_height=float(params.get("figureHeight", 6.0)),
        dpi=int(params.get("dpi", 300)),
    )

    return {"graph": graphPath, "data": dataPath}


VOLCANO_BLOCK = BlockDefinition(
    key="volcanoPlot",
    label="Volcano Plot",
    configSection="volcanoPlot",
    inputs=lambda node: [
        PortSpec("file1", "File 1", C.UNIFIER),
        PortSpec("file2", "File 2", C.UNIFIER),
        PortSpec("ratio", "Ratio", C.LITERAL),
        PortSpec("pvalue", "P-Value", C.LITERAL),
    ],
    outputs=lambda node: [
        PortSpec("graph", "Graph Plot", C.OUTPUT_PNG),
        PortSpec("data", "Data", C.UNIFIER),
    ],
    execute=_executeVolcano,
    description="Compares two unified runs by ratio and p-value, exporting one quadrant.",
)


# --------------------------------------------------------------------------- #
#  Upset Plot                                                                  #
# --------------------------------------------------------------------------- #

def _executeUpset(node, inputs, ctx):
    from CLI_upsetplot import run_upset_plot

    paths = _collectDynamic(node, inputs, "Upset Plot")

    params = dict(node.params)

    #An empty selection means every file takes part, which is the intersection
    #across all of them.
    selection = params.get("selection") or []
    selection = [bool(flag) for flag in selection][:len(paths)]
    selection = selection + [True] * (len(paths) - len(selection))

    graphPath = ctx.path(node, "upset.png")
    dataPath = ctx.path(node, "upsetIntersection.csv")

    run_upset_plot(
        files=paths,
        graph_output=graphPath,
        export_intersection=True,
        intersection_output=dataPath,
        selection=selection,
        figure_width=float(params.get("figureWidth", 14.0)),
        figure_height=float(params.get("figureHeight", 6.0)),
        dpi=int(params.get("dpi", 300)),
    )

    return {"plot": graphPath, "data": dataPath}


UPSET_BLOCK = BlockDefinition(
    key="upsetPlot",
    label="Upset Plot",
    configSection="upsetPlot",
    #inputCount is driven by the + and - buttons on the block, not the popup
    configKeys=["selection", "figureWidth", "figureHeight", "dpi"],
    inputs=lambda node: _dynamicInputs(node, C.UNIFIER, "File"),
    outputs=lambda node: [
        PortSpec("plot", "Plot", C.OUTPUT_PNG),
        PortSpec("data", "Data", C.UNIFIER),
    ],
    execute=_executeUpset,
    dynamicInputs=True,
    description="Shows how the sequences of several runs overlap and exports the "
                "selected intersection.",
)


# --------------------------------------------------------------------------- #
#  Ranking Plot                                                                #
# --------------------------------------------------------------------------- #

def _executeRanking(node, inputs, ctx):
    from CLI_rankingPlot import run_ranking_plot

    fileA = _require(inputs, "fileA", "Ranking Plot", "File A")
    fileB = _require(inputs, "fileB", "Ranking Plot", "File B")

    params = dict(node.params)

    #Unconnected target ports are skipped; the plot only boxes what is wired in
    targets = [str(inputs[spec.key].value) for spec in node.inputs()
               if spec.key.startswith("target") and inputs.get(spec.key) is not None]

    topN = value(inputs, "topN", params.get("points", 100))
    slope = value(inputs, "slope", params.get("slope", 1.0))

    graphPath = ctx.path(node, "ranking.png")
    rankingPath = ctx.path(node, "ranking.csv")
    abovePath = ctx.path(node, "aboveLine.csv")

    run_ranking_plot(
        fileA, fileB,
        output_path=graphPath,
        slope=float(slope),
        b=float(params.get("b", 0.0)),
        points=float(topN),
        percent_or_count=params.get("percentOrCount", "#"),
        count_file1=bool(params.get("countFile1", True)),
        count_file2=bool(params.get("countFile2", False)),
        export_graph=True,
        log_scale=bool(params.get("logScale", False)),
        ranking_output=rankingPath,
        above_output=abovePath,
        figure_width=float(params.get("figureWidth", 6.0)),
        figure_height=float(params.get("figureHeight", 4.0)),
        dpi=int(params.get("dpi", 300)),
        target_sequences=targets,
    )

    return {"above": abovePath, "graph": graphPath, "ranking": rankingPath}


RANKING_BLOCK = BlockDefinition(
    key="rankingPlot",
    label="Ranking Plot",
    configSection="rankingPlot",
    configKeys=["slope", "b", "points", "percentOrCount", "countFile1",
                "countFile2", "logScale", "figureWidth", "figureHeight", "dpi"],
    inputs=lambda node: [
        PortSpec("fileA", "File A", C.UNIFIER),
        PortSpec("fileB", "File B", C.UNIFIER),
        PortSpec("topN", "Top N", C.LITERAL),
        PortSpec("slope", "Slope (M)", C.LITERAL),
        *[PortSpec(f"target{i + 1}", f"Target {i + 1}", C.LITERAL)
          for i in range(max(1, node.inputCount))],
    ],
    outputs=lambda node: [
        PortSpec("above", "Above-Line Data", C.UNIFIER),
        PortSpec("graph", "Graph Plot", C.OUTPUT_PNG),
        PortSpec("ranking", "Ranking", C.RANKING_OUTPUT),
    ],
    execute=_executeRanking,
    dynamicInputs=True,
    dynamicLabel="Target",
    description="Compares two unified runs by rank, splitting them on a line.",
)


# --------------------------------------------------------------------------- #
#  Registry                                                                    #
# --------------------------------------------------------------------------- #

REGISTRY = {
    block.key: block
    for block in (
        FILE_INPUT_BLOCK,
        LITERAL_BLOCK,
        FASTA_BLOCK,
        PAIRED_ASSEMBLER_BLOCK,
        SEQUENCE_COUNTER_BLOCK,
        RUN_UNIFIER_BLOCK,
        VOLCANO_BLOCK,
        UPSET_BLOCK,
        RANKING_BLOCK,
        OUTPUT_BLOCK,
    )
}

#The order blocks appear in the palette, grouped from inputs through to outputs
PALETTE_ORDER = (
    ("Inputs", ("fileInput", "literal")),
    ("Ingestion", ("fasta", "pairedAssembler")),
    ("Counting", ("sequenceCounter", "runUnifier")),
    ("Plots", ("volcanoPlot", "upsetPlot", "rankingPlot")),
    ("Outputs", ("output",)),
)
