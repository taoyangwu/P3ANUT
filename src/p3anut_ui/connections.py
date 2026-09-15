"""
Connection types for the P3ANUT node graph.

A connection's shape encodes how far through the pipeline its data has
travelled: a raw File Input is a circle, the output of the Paired Assembler and
the FASTA block is a triangle, and the side count climbs from there. A two
sided shape is skipped so nothing is ever mistaken for a plain line.

Colours come from the matplotlib "tab" (Tableau) cycle, so a connection is
identifiable by colour alone as well as by shape.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class ConnectionType:
    key: str
    label: str
    sides: int          # 0 draws a circle; 3 and up draw a regular polygon
    color: str          # matplotlib tab cycle, as a hex value
    extension: str      # file extension an Output block writes for this type
    stage: int          # pipeline stages the data has passed through


#--------------------------------------------------------------------------#
# Stage 0 - raw inputs the user supplies directly
#--------------------------------------------------------------------------#
FILE_INPUT = ConnectionType("file", "File Input", 0, "#1f77b4", "", 0)
LITERAL = ConnectionType("literal", "Literal Input", 0, "#ff7f0e", ".txt", 0)

#--------------------------------------------------------------------------#
# Stage 1 - ingestion, the Paired Assembler and the FASTA block
#--------------------------------------------------------------------------#
JSON = ConnectionType("json", "JSON", 3, "#2ca02c", ".json", 1)
METADATA = ConnectionType("metadata", "Meta Data", 3, "#d62728", ".json", 1)

#--------------------------------------------------------------------------#
# Stage 2 and beyond - counting, unifying and the plots
#--------------------------------------------------------------------------#
SEQUENCE_COUNT = ConnectionType("sequenceCount", "Sequence Counter", 4, "#9467bd", ".csv", 2)
UNIFIER = ConnectionType("unifier", "Unifier", 5, "#8c564b", ".csv", 3)
OUTPUT_PNG = ConnectionType("png", "Output PNG", 6, "#e377c2", ".png", 4)
RANKING_OUTPUT = ConnectionType("ranking", "Ranking Output", 7, "#7f7f7f", ".csv", 5)

#--------------------------------------------------------------------------#
# The Output block accepts and forwards anything, so its ports carry no
# stage of their own and take on whatever is wired into them.
#--------------------------------------------------------------------------#
ANY = ConnectionType("any", "Any", 0, "#bcbd22", "", 0)


ALL_TYPES = (FILE_INPUT, LITERAL, JSON, METADATA, SEQUENCE_COUNT,
             UNIFIER, OUTPUT_PNG, RANKING_OUTPUT, ANY)

BY_KEY = {t.key: t for t in ALL_TYPES}


def get(key):
    """Look a connection type up by its key."""
    return BY_KEY[key]


def compatible(sourceType, targetType):
    """
    Whether an edge may run from an output of sourceType into an input of
    targetType. ANY, which only the Output block uses, pairs with everything.
    """
    if ANY in (sourceType, targetType):
        return True

    return sourceType.key == targetType.key


def extensionFor(connectionType, sourcePath=None):
    """
    The file extension an Output block should use for this connection type.

    A raw File Input has no extension of its own, so it keeps the one the
    selected file already had.
    """
    if connectionType.extension:
        return connectionType.extension

    if sourcePath:
        import os
        return os.path.splitext(sourcePath)[1] or ".dat"

    return ".dat"
