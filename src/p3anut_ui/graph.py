"""
The node graph model: blocks, the connections between them, and the
user-readable file the whole thing is saved to.

Saving exists purely so a pipeline does not have to be rebuilt by hand.
Reopening a saved file restores every block, connection and configured value so
the user can keep editing and run it again. It says nothing about a previous
run: no output is cached and no partially finished run is resumed.
"""

import itertools
import os

import yaml

from . import blocks as B
from . import connections as C
from .config import ConfigTemplate

FILE_VERSION = 1


class Node:
    """One block placed on the canvas, holding its own copy of the parameters."""

    _counter = itertools.count(1)

    def __init__(self, blockKey, x=0.0, y=0.0, nodeId=None, params=None,
                 inputCount=None, template=None):
        self.blockKey = blockKey
        self.definition = B.REGISTRY[blockKey]
        self.id = nodeId or f"{blockKey}_{next(Node._counter)}"
        self.x = float(x)
        self.y = float(y)

        self.inputCount = (self._defaultInputCount(template)
                           if inputCount is None else int(inputCount))

        #Each instance gets its own copy of the template defaults, so editing
        #this block's popup never touches another block of the same type.
        self.params = dict(params) if params is not None else self._seedParams(template)

    def _defaultInputCount(self, template):
        if not self.definition.dynamicInputs:
            return 0

        if template is not None:
            declared = template.parameter(self.definition.configSection, "inputCount")
            if declared and declared.get("Value"):
                return int(declared["Value"])

        return 2

    def _seedParams(self, template):
        if template is None:
            return {}

        return template.defaults(self.definition.configSection,
                                 self.definition.configKeys)

    @property
    def label(self):
        return self.definition.label

    def inputs(self):
        return self.definition.inputs(self)

    def outputs(self):
        return self.definition.outputs(self)

    def port(self, key, isInput):
        for spec in (self.inputs() if isInput else self.outputs()):
            if spec.key == key:
                return spec

        return None


class Edge:
    """A connection from one block's output to another block's input."""

    def __init__(self, sourceNode, sourcePort, targetNode, targetPort):
        self.sourceNode = sourceNode
        self.sourcePort = sourcePort
        self.targetNode = targetNode
        self.targetPort = targetPort

    def key(self):
        return (self.sourceNode, self.sourcePort, self.targetNode, self.targetPort)

    def __eq__(self, other):
        return isinstance(other, Edge) and self.key() == other.key()

    def __hash__(self):
        return hash(self.key())


class Graph:
    """The blocks, their connections and the rules that govern them."""

    def __init__(self, template=None):
        self.template = template or ConfigTemplate()
        self.nodes = {}
        self.edges = []

    # ---------------------------------------------------------------- blocks
    def addNode(self, blockKey, x=0.0, y=0.0, **kwargs):
        node = Node(blockKey, x, y, template=self.template, **kwargs)
        self.nodes[node.id] = node
        return node

    def removeNode(self, nodeId):
        self.nodes.pop(nodeId, None)
        self.edges = [e for e in self.edges
                      if e.sourceNode != nodeId and e.targetNode != nodeId]

    def setInputCount(self, nodeId, count):
        """Grow or shrink a dynamic block, dropping edges to removed ports."""
        node = self.nodes[nodeId]
        count = max(1, int(count))
        node.inputCount = count

        valid = {spec.key for spec in node.inputs()}
        self.edges = [e for e in self.edges
                      if e.targetNode != nodeId or e.targetPort in valid]

    # ----------------------------------------------------------------- edges
    def incoming(self, nodeId, portKey=None):
        return [e for e in self.edges
                if e.targetNode == nodeId and (portKey is None or e.targetPort == portKey)]

    def outgoing(self, nodeId, portKey=None):
        return [e for e in self.edges
                if e.sourceNode == nodeId and (portKey is None or e.sourcePort == portKey)]

    def resolvedOutputType(self, nodeId, portKey, _seen=None):
        """
        The connection type an output actually carries.

        The Output block is a pass-through whose ports are declared as ANY, so
        its real type is whatever is wired into it. Everything else simply
        reports the type its port declares.
        """
        node = self.nodes[nodeId]
        spec = node.port(portKey, isInput=False)

        if spec is None:
            return C.ANY

        if spec.type is not C.ANY:
            return spec.type

        _seen = _seen or set()
        if nodeId in _seen:
            return C.ANY
        _seen.add(nodeId)

        feeding = self.incoming(nodeId, "data")
        if not feeding:
            return C.ANY

        edge = feeding[0]
        return self.resolvedOutputType(edge.sourceNode, edge.sourcePort, _seen)

    def canConnect(self, sourceNode, sourcePort, targetNode, targetPort):
        """
        Whether an edge is allowed, and if not, why. The reason is shown to the
        user rather than the connection just silently refusing to form.
        """
        if sourceNode == targetNode:
            return False, "A block cannot be connected to itself"

        source = self.nodes.get(sourceNode)
        target = self.nodes.get(targetNode)

        if source is None or target is None:
            return False, "Unknown block"

        sourceSpec = source.port(sourcePort, isInput=False)
        targetSpec = target.port(targetPort, isInput=True)

        if sourceSpec is None or targetSpec is None:
            return False, "Connections run from an output on the right to an input on the left"

        #An input takes a single connection; an output may feed many
        if self.incoming(targetNode, targetPort):
            return False, f"'{targetSpec.label}' already has a connection"

        sourceType = self.resolvedOutputType(sourceNode, sourcePort)

        #A Ranking Output is a terminal record of one step - it only ever goes
        #into an Output block.
        if sourceType is C.RANKING_OUTPUT and target.blockKey != "output":
            return False, "A Ranking Output can only be connected to an Output block"

        if not C.compatible(sourceType, targetSpec.type):
            return False, (f"{sourceType.label} cannot connect to "
                           f"{targetSpec.type.label}")

        if self._wouldCycle(sourceNode, targetNode):
            return False, "That connection would create a loop"

        return True, ""

    def _wouldCycle(self, sourceNode, targetNode):
        """True when targetNode already feeds sourceNode, directly or not."""
        stack = [targetNode]
        seen = set()

        while stack:
            current = stack.pop()
            if current == sourceNode:
                return True
            if current in seen:
                continue
            seen.add(current)
            stack.extend(e.targetNode for e in self.outgoing(current))

        return False

    def addEdge(self, sourceNode, sourcePort, targetNode, targetPort):
        allowed, reason = self.canConnect(sourceNode, sourcePort, targetNode, targetPort)
        if not allowed:
            raise ValueError(reason)

        edge = Edge(sourceNode, sourcePort, targetNode, targetPort)
        self.edges.append(edge)
        return edge

    def removeEdge(self, edge):
        if edge in self.edges:
            self.edges.remove(edge)

    # ------------------------------------------------------------- execution
    def executionOrder(self):
        """
        The blocks in dependency order.

        Raises when the graph contains a loop, which the connection rules
        already prevent, but a hand-edited save file could still carry one.
        """
        pending = {nodeId: {e.sourceNode for e in self.incoming(nodeId)}
                   for nodeId in self.nodes}

        ordered = []
        while pending:
            ready = sorted(nodeId for nodeId, deps in pending.items() if not deps)

            if not ready:
                raise ValueError("The graph contains a loop and cannot be run")

            for nodeId in ready:
                ordered.append(nodeId)
                pending.pop(nodeId)

            for deps in pending.values():
                deps.difference_update(ready)

        return ordered

    # ------------------------------------------------------------ save, load
    def toDocument(self):
        return {
            "version": FILE_VERSION,
            "nodes": [
                {
                    "id": node.id,
                    "block": node.blockKey,
                    "x": round(node.x, 2),
                    "y": round(node.y, 2),
                    **({"inputCount": node.inputCount}
                       if node.definition.dynamicInputs else {}),
                    "params": node.params,
                }
                for node in self.nodes.values()
            ],
            "edges": [
                {
                    "from": {"block": e.sourceNode, "port": e.sourcePort},
                    "to": {"block": e.targetNode, "port": e.targetPort},
                }
                for e in self.edges
            ],
        }

    def save(self, path):
        document = self.toDocument()

        with open(path, "w") as handle:
            handle.write("# P3ANUT pipeline graph\n"
                         "# Blocks, their connections and their configured values.\n"
                         "# Reopen this file to keep editing the pipeline and run it again.\n\n")
            yaml.safe_dump(document, handle, sort_keys=False, default_flow_style=False)

        return path

    @classmethod
    def load(cls, path, template=None):
        with open(path, "r") as handle:
            document = yaml.safe_load(handle) or {}

        version = document.get("version")
        if version != FILE_VERSION:
            raise ValueError(f"Unsupported graph file version: {version}")

        graph = cls(template)

        for record in document.get("nodes", []):
            blockKey = record.get("block")
            if blockKey not in B.REGISTRY:
                raise ValueError(f"Unknown block type in saved graph: {blockKey}")

            node = Node(
                blockKey,
                record.get("x", 0.0),
                record.get("y", 0.0),
                nodeId=record.get("id"),
                params=record.get("params") or {},
                inputCount=record.get("inputCount"),
                template=graph.template,
            )
            graph.nodes[node.id] = node

        for record in document.get("edges", []):
            source = record.get("from", {})
            target = record.get("to", {})
            graph.edges.append(Edge(source.get("block"), source.get("port"),
                                    target.get("block"), target.get("port")))

        return graph
