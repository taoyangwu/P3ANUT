"""
The single P3ANUT window: the canvas on the left, the block palette on the
right, and the Run button at the top of that palette.
"""

import os

from PyQt6.QtCore import QMimeData, QSize, Qt, pyqtSignal
from PyQt6.QtGui import QAction, QColor, QDrag, QPainter, QPixmap
from PyQt6.QtWidgets import (
    QFileDialog,
    QFrame,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QPlainTextEdit,
    QPushButton,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from . import blocks as B
from . import connections as C
from .canvas import BLOCK_MIME, Canvas
from .config import ConfigTemplate
from .executor import GraphRunner, RunThread
from .graph import Graph
from .items import portShapePath

GRAPH_FILTER = "P3ANUT pipeline (*.p3g.yaml *.yaml);;All files (*)"


class BlockPalette(QListWidget):
    """The draggable list of available blocks."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setDragEnabled(True)
        self.setIconSize(QSize(14, 14))
        self.setAlternatingRowColors(False)
        self.setStyleSheet("QListWidget { border: none; } "
                           "QListWidget::item { padding: 5px 4px; }")

        for group, keys in B.PALETTE_ORDER:
            header = QListWidgetItem(group.upper())
            header.setFlags(Qt.ItemFlag.NoItemFlags)
            header.setForeground(QColor("#888888"))
            self.addItem(header)

            for key in keys:
                definition = B.REGISTRY[key]
                item = QListWidgetItem(f"  {definition.label}")
                item.setData(Qt.ItemDataRole.UserRole, key)
                item.setToolTip(definition.description)
                self.addItem(item)

    def startDrag(self, actions):
        item = self.currentItem()
        blockKey = item.data(Qt.ItemDataRole.UserRole) if item else None

        if not blockKey:
            return

        mime = QMimeData()
        mime.setData(BLOCK_MIME, blockKey.encode())

        drag = QDrag(self)
        drag.setMimeData(mime)
        drag.exec(Qt.DropAction.CopyAction)


class ConnectionLegend(QWidget):
    """Shows what each connection shape and colour means."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumHeight(len(C.ALL_TYPES) * 18 + 8)

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)

        for index, connectionType in enumerate(C.ALL_TYPES):
            y = 12 + index * 18

            painter.save()
            painter.translate(12, y)
            painter.setBrush(QColor(connectionType.color))
            painter.setPen(QColor("#33000000"))
            painter.drawPath(portShapePath(connectionType.sides, 6))
            painter.restore()

            painter.setPen(QColor("#444444"))
            shape = "circle" if connectionType.sides < 3 else f"{connectionType.sides} sides"
            painter.drawText(28, y + 4, f"{connectionType.label}  ({shape})")


class MainWindow(QMainWindow):
    """The whole application."""

    def __init__(self, graphPath=None):
        super().__init__()

        self.template = ConfigTemplate()
        self.graph = Graph(self.template)
        self.currentPath = None
        self.runThread = None
        self.runner = None

        self.setWindowTitle("P3ANUT")
        self.resize(1440, 900)

        self.canvas = Canvas(self.graph, self)
        self.canvas.statusMessage.connect(self.showStatus)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(self.canvas)
        splitter.addWidget(self._buildSidebar())
        splitter.setStretchFactor(0, 1)
        splitter.setSizes([1100, 340])

        self.setCentralWidget(splitter)
        self._buildMenu()

        self.statusBar().showMessage("Drag a block from the right onto the canvas to begin")

        if graphPath:
            self.loadGraph(graphPath)
        else:
            self.canvas.rebuild()

    # -------------------------------------------------------------- sidebar
    def _buildSidebar(self):
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)

        #Run controls sit at the very top of the right hand bar
        self.runButton = QPushButton("▶  Run")
        self.runButton.setMinimumHeight(38)
        self.runButton.setStyleSheet(
            "QPushButton { background: #2ca02c; color: white; font-weight: bold;"
            " border: none; border-radius: 5px; font-size: 14px; }"
            "QPushButton:hover { background: #249024; }"
            "QPushButton:disabled { background: #a8c8a8; }")
        self.runButton.clicked.connect(self.runGraph)
        layout.addWidget(self.runButton)

        self.progress = QProgressBar()
        self.progress.setMinimum(0)
        self.progress.setMaximum(100)
        self.progress.setValue(0)
        self.progress.setFormat("idle")
        layout.addWidget(self.progress)

        divider = QFrame()
        divider.setFrameShape(QFrame.Shape.HLine)
        divider.setStyleSheet("color: #dddddd;")
        layout.addWidget(divider)

        label = QLabel("Blocks")
        label.setStyleSheet("font-weight: bold;")
        layout.addWidget(label)

        self.palette = BlockPalette()
        layout.addWidget(self.palette, 2)

        legendLabel = QLabel("Connections")
        legendLabel.setStyleSheet("font-weight: bold;")
        layout.addWidget(legendLabel)
        layout.addWidget(ConnectionLegend())

        logLabel = QLabel("Run log")
        logLabel.setStyleSheet("font-weight: bold;")
        layout.addWidget(logLabel)

        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(500)
        self.log.setStyleSheet("font-family: Menlo, Consolas, 'Courier New', monospace;"
                               " font-size: 11px;")
        layout.addWidget(self.log, 1)

        return panel

    def _buildMenu(self):
        fileMenu = self.menuBar().addMenu("&File")

        newAction = QAction("&New pipeline", self)
        newAction.setShortcut("Ctrl+N")
        newAction.triggered.connect(self.newGraph)
        fileMenu.addAction(newAction)

        openAction = QAction("&Open pipeline...", self)
        openAction.setShortcut("Ctrl+O")
        openAction.triggered.connect(self.openGraph)
        fileMenu.addAction(openAction)

        saveAction = QAction("&Save pipeline", self)
        saveAction.setShortcut("Ctrl+S")
        saveAction.triggered.connect(self.saveGraph)
        fileMenu.addAction(saveAction)

        saveAsAction = QAction("Save pipeline &as...", self)
        saveAsAction.setShortcut("Ctrl+Shift+S")
        saveAsAction.triggered.connect(lambda: self.saveGraph(forcePrompt=True))
        fileMenu.addAction(saveAsAction)

        fileMenu.addSeparator()

        quitAction = QAction("&Quit", self)
        quitAction.setShortcut("Ctrl+Q")
        quitAction.triggered.connect(self.close)
        fileMenu.addAction(quitAction)

        editMenu = self.menuBar().addMenu("&Edit")
        deleteAction = QAction("&Delete selection", self)
        deleteAction.setShortcut("Delete")
        deleteAction.triggered.connect(self.canvas.deleteSelection)
        editMenu.addAction(deleteAction)

    # ------------------------------------------------------------ file menu
    def newGraph(self):
        self.graph = Graph(self.template)
        self.canvas.graph = self.graph
        self.currentPath = None
        self.canvas.rebuild()
        self.setWindowTitle("P3ANUT")
        self.showStatus("New pipeline")

    def openGraph(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open a pipeline", os.path.expanduser("~"), GRAPH_FILTER)

        if path:
            self.loadGraph(path)

    def loadGraph(self, path):
        try:
            self.graph = Graph.load(path, self.template)
        except Exception as exc:                          # noqa: BLE001 - shown to the user
            QMessageBox.critical(self, "Could not open pipeline", str(exc))
            return

        self.canvas.graph = self.graph
        self.currentPath = path
        self.canvas.rebuild()
        self.setWindowTitle(f"P3ANUT - {os.path.basename(path)}")
        self.showStatus(f"Opened {path}")

    def saveGraph(self, forcePrompt=False):
        path = self.currentPath

        if forcePrompt or not path:
            path, _ = QFileDialog.getSaveFileName(
                self, "Save the pipeline",
                os.path.join(os.path.expanduser("~"), "pipeline.p3g.yaml"),
                GRAPH_FILTER)

            if not path:
                return

        try:
            self.graph.save(path)
        except Exception as exc:                          # noqa: BLE001 - shown to the user
            QMessageBox.critical(self, "Could not save pipeline", str(exc))
            return

        self.currentPath = path
        self.setWindowTitle(f"P3ANUT - {os.path.basename(path)}")
        self.showStatus(f"Saved to {path}")

    # ---------------------------------------------------------- running
    def runGraph(self):
        if self.runThread is not None:
            return

        if not self.graph.nodes:
            QMessageBox.information(self, "Nothing to run",
                                    "Drag some blocks onto the canvas first.")
            return

        self.log.clear()
        self.runButton.setEnabled(False)
        self.progress.setValue(0)
        self.progress.setFormat("starting...")

        self.runner = GraphRunner(self.graph)
        self.runner.progress.connect(self.onProgress)
        self.runner.message.connect(self.onMessage)
        self.runner.finished.connect(self.onFinished)

        self.runThread = RunThread(self.runner, self)
        self.runThread.start()

    def onProgress(self, completed, total, label):
        #Progress is completed blocks over total blocks in the graph
        self.progress.setMaximum(total)
        self.progress.setValue(completed)
        self.progress.setFormat(f"{completed}/{total} - {label}")

    def onMessage(self, message):
        self.log.appendPlainText(message)

    def onFinished(self, succeeded, summary):
        self.runButton.setEnabled(True)
        self.log.appendPlainText(summary)

        if self.runThread is not None:
            self.runThread.quit()
            self.runThread.wait()
            self.runThread = None
        self.runner = None

        if succeeded:
            self.progress.setFormat("complete")
            self.showStatus(summary)
        else:
            self.progress.setFormat("failed")
            QMessageBox.critical(self, "The run stopped", summary)

    # ----------------------------------------------------------------- misc
    def showStatus(self, message):
        self.statusBar().showMessage(message, 6000)

    def closeEvent(self, event):
        if self.runThread is not None:
            self.runner.cancel()
            self.runThread.quit()
            self.runThread.wait(3000)

        super().closeEvent(event)
