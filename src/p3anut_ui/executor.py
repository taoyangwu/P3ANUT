"""
Running a graph.

Every block in the graph is executed in dependency order. Progress is reported
as completed blocks over total blocks. If any block raises, the run stops there
and then - unaffected branches are not carried on with, and nothing is resumed
on the next run.
"""

import os
import shutil
import tempfile
import traceback

import matplotlib
#Figures are rendered off the GUI thread, so the interactive backend must not
#be involved.
matplotlib.use("Agg", force=True)

from PyQt6.QtCore import QObject, QThread, pyqtSignal

from . import blocks as B


class RunContext:
    """Per-run scratch space and logging handed to each block's execute."""

    def __init__(self, workdir, logger=None):
        self.workdir = workdir
        self._logger = logger

    def path(self, node, name):
        """A path inside this run's working directory, unique to the block."""
        directory = os.path.join(self.workdir, node.id)
        os.makedirs(directory, exist_ok=True)
        return os.path.join(directory, name)

    def log(self, message):
        if self._logger:
            self._logger(message)
        else:
            print(message)


class GraphRunner(QObject):
    """Executes a graph, reporting progress and the first failure."""

    progress = pyqtSignal(int, int, str)     # completed, total, current block label
    message = pyqtSignal(str)
    finished = pyqtSignal(bool, str)         # succeeded, summary

    def __init__(self, graph, keepWorkdir=False):
        super().__init__()
        self.graph = graph
        self.keepWorkdir = keepWorkdir
        self.workdir = None
        self._cancelled = False

    def cancel(self):
        self._cancelled = True

    def run(self):
        total = len(self.graph.nodes)

        if total == 0:
            self.finished.emit(False, "There is nothing on the canvas to run.")
            return

        self.workdir = tempfile.mkdtemp(prefix="p3anut_run_")
        context = RunContext(self.workdir, self.message.emit)

        self.message.emit(f"Run started - working directory {self.workdir}")

        #Outputs already produced, keyed by (blockId, portKey)
        produced = {}
        completed = 0

        try:
            order = self.graph.executionOrder()
        except ValueError as exc:
            self.finished.emit(False, str(exc))
            return

        for nodeId in order:
            if self._cancelled:
                self.finished.emit(False, "Run cancelled.")
                return

            node = self.graph.nodes[nodeId]
            self.progress.emit(completed, total, node.label)

            inputs = {}
            for edge in self.graph.incoming(nodeId):
                key = (edge.sourceNode, edge.sourcePort)
                if key in produced:
                    inputs[edge.targetPort] = produced[key]

            try:
                results = node.definition.execute(node, inputs, context) or {}
            except Exception as exc:                      # noqa: BLE001 - reported to the user
                self.message.emit(traceback.format_exc())
                self._cleanup(failed=True)
                self.finished.emit(
                    False,
                    f"{node.label} failed: {exc}\n\nThe run stopped at that block.")
                return

            for spec in node.outputs():
                if spec.key in results:
                    produced[(nodeId, spec.key)] = B.Payload(
                        results[spec.key],
                        self.graph.resolvedOutputType(nodeId, spec.key),
                    )

            completed += 1
            self.progress.emit(completed, total, node.label)

        self._cleanup(failed=False)
        self.finished.emit(True, f"Run complete - {completed} of {total} blocks executed.")

    def _cleanup(self, failed):
        """
        Remove the run's scratch directory.

        It is kept after a failure, and whenever the user asked for it, because
        the partial results are the most useful thing to look at when a block
        has gone wrong.
        """
        if self.keepWorkdir or failed:
            self.message.emit(f"Working directory kept at {self.workdir}")
            return

        shutil.rmtree(self.workdir, ignore_errors=True)


class RunThread(QThread):
    """Carries a GraphRunner off the GUI thread so the window stays responsive."""

    def __init__(self, runner, parent=None):
        super().__init__(parent)
        self.runner = runner
        self.runner.moveToThread(self)

    def run(self):
        self.runner.run()
