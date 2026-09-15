"""
The canvas: the left hand area where blocks are placed and wired together.

Blocks are dragged in from the palette on the right. Dragging from one port to
another makes a connection, which is only formed when the two ports are
compatible - otherwise the reason is reported rather than the connection
quietly failing.
"""

from PyQt6.QtCore import QPointF, Qt, pyqtSignal
from PyQt6.QtGui import QColor, QPainter, QPainterPath, QPen
from PyQt6.QtWidgets import QGraphicsPathItem, QGraphicsScene, QGraphicsView

from .dialogs import ParameterDialog
from .items import EdgeItem, NodeItem, PortItem

BLOCK_MIME = "application/x-p3anut-block"


class Canvas(QGraphicsView):
    """The node editing surface."""

    statusMessage = pyqtSignal(str)
    graphChanged = pyqtSignal()

    def __init__(self, graph, parent=None):
        super().__init__(parent)
        self.graph = graph
        self.nodeItems = {}
        self.edgeItems = []

        self._pendingPort = None
        self._tempEdge = None

        scene = QGraphicsScene(self)
        scene.setSceneRect(-2000, -2000, 6000, 6000)
        self.setScene(scene)

        self.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        self.setDragMode(QGraphicsView.DragMode.RubberBandDrag)
        self.setAcceptDrops(True)
        self.setBackgroundBrush(QColor("#f2f3f5"))
        self.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)

    # ------------------------------------------------------------- building
    def rebuild(self):
        """Redraw the whole canvas from the graph model."""
        self.scene().clear()
        self.nodeItems.clear()
        self.edgeItems.clear()
        self._pendingPort = None
        self._tempEdge = None

        for node in self.graph.nodes.values():
            item = NodeItem(self, node)
            self.scene().addItem(item)
            self.nodeItems[node.id] = item

        for edge in list(self.graph.edges):
            sourceItem = self.nodeItems.get(edge.sourceNode)
            targetItem = self.nodeItems.get(edge.targetNode)

            if sourceItem is None or targetItem is None:
                self.graph.removeEdge(edge)
                continue

            sourcePort = self._findPort(sourceItem, edge.sourcePort, isInput=False)
            targetPort = self._findPort(targetItem, edge.targetPort, isInput=True)

            if sourcePort is None or targetPort is None:
                self.graph.removeEdge(edge)
                continue

            item = EdgeItem(self, edge, sourcePort, targetPort)
            self.scene().addItem(item)
            self.edgeItems.append(item)

        self.refreshTypes()
        self.graphChanged.emit()

    @staticmethod
    def _findPort(nodeItem, key, isInput):
        ports = nodeItem.inputPorts if isInput else nodeItem.outputPorts
        for port in ports:
            if port.spec.key == key:
                return port
        return None

    def addBlock(self, blockKey, scenePos):
        node = self.graph.addNode(blockKey, scenePos.x(), scenePos.y())
        item = NodeItem(self, node)
        self.scene().addItem(item)
        self.nodeItems[node.id] = item

        self.statusMessage.emit(f"Added {node.label}")
        self.graphChanged.emit()
        return item

    # --------------------------------------------------------------- typing
    def refreshTypes(self):
        """Re-resolve pass-through port types and recolour everything."""
        for item in self.nodeItems.values():
            item.refreshPortTypes()

        for edge in self.edgeItems:
            edge.refresh()

    def refreshEdges(self):
        for edge in self.edgeItems:
            edge.refresh()

    # ---------------------------------------------------------- connections
    def beginConnection(self, port):
        """Start dragging a new connection out of a port."""
        self._pendingPort = port

        self._tempEdge = QGraphicsPathItem()
        self._tempEdge.setPen(QPen(QColor(port.connectionType().color), 2,
                                   Qt.PenStyle.DashLine))
        self._tempEdge.setZValue(5)
        self.scene().addItem(self._tempEdge)

    def _updateTempEdge(self, scenePos):
        start = self._pendingPort.scenePos()
        reach = max(40.0, abs(scenePos.x() - start.x()) * 0.5)

        path = QPainterPath(start)
        if self._pendingPort.isInput:
            path.cubicTo(QPointF(start.x() - reach, start.y()),
                         QPointF(scenePos.x() + reach, scenePos.y()), scenePos)
        else:
            path.cubicTo(QPointF(start.x() + reach, start.y()),
                         QPointF(scenePos.x() - reach, scenePos.y()), scenePos)

        self._tempEdge.setPath(path)

    def _finishConnection(self, scenePos):
        source = self._pendingPort
        target = self._portAt(scenePos)

        self.scene().removeItem(self._tempEdge)
        self._tempEdge = None
        self._pendingPort = None

        if target is None or target is source:
            return

        #A connection always runs output to input, whichever end was grabbed first
        if source.isInput and not target.isInput:
            source, target = target, source

        if source.isInput or not target.isInput:
            self.statusMessage.emit("Connect an output on the right to an input on the left")
            return

        try:
            edge = self.graph.addEdge(source.nodeItem.node.id, source.spec.key,
                                      target.nodeItem.node.id, target.spec.key)
        except ValueError as exc:
            self.statusMessage.emit(str(exc))
            return

        item = EdgeItem(self, edge, source, target)
        self.scene().addItem(item)
        self.edgeItems.append(item)

        self.refreshTypes()
        self.statusMessage.emit(
            f"Connected {source.spec.label} to {target.spec.label}")
        self.graphChanged.emit()

    def _portAt(self, scenePos):
        for item in self.scene().items(scenePos):
            if isinstance(item, PortItem):
                return item
        return None

    # ------------------------------------------------------------- dynamics
    def changeInputCount(self, nodeItem, delta):
        node = nodeItem.node
        newCount = max(1, node.inputCount + delta)

        if newCount == node.inputCount:
            self.statusMessage.emit(f"{node.label} needs at least one input")
            return

        self.graph.setInputCount(node.id, newCount)

        #Edges to removed ports are gone from the model; redraw to match
        self.rebuild()
        self.statusMessage.emit(f"{node.label} now has {newCount} file inputs")

    def openParameters(self, nodeItem):
        dialog = ParameterDialog(nodeItem.node, self.graph.template, self)

        if dialog.exec():
            dialog.applyTo(nodeItem.node)
            #A changed subtitle moves the port rows, so redo the layout
            nodeItem.refreshLayout()
            self.refreshEdges()
            self.statusMessage.emit(f"Updated {nodeItem.node.label} parameters")
            self.graphChanged.emit()

    def deleteSelection(self):
        removed = 0

        for item in list(self.scene().selectedItems()):
            if isinstance(item, EdgeItem):
                self.graph.removeEdge(item.edge)
                removed += 1
            elif isinstance(item, NodeItem):
                self.graph.removeNode(item.node.id)
                removed += 1

        if removed:
            self.rebuild()
            self.statusMessage.emit(f"Removed {removed} item(s)")

    # ---------------------------------------------------------------- input
    def mouseMoveEvent(self, event):
        if self._tempEdge is not None:
            self._updateTempEdge(self.mapToScene(event.pos()))
            event.accept()
            return

        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if self._tempEdge is not None:
            self._finishConnection(self.mapToScene(event.pos()))
            event.accept()
            return

        super().mouseReleaseEvent(event)

    def keyPressEvent(self, event):
        if event.key() in (Qt.Key.Key_Delete, Qt.Key.Key_Backspace):
            self.deleteSelection()
            event.accept()
            return

        super().keyPressEvent(event)

    def wheelEvent(self, event):
        if event.modifiers() & Qt.KeyboardModifier.ControlModifier:
            factor = 1.15 if event.angleDelta().y() > 0 else 1 / 1.15
            self.scale(factor, factor)
            event.accept()
            return

        super().wheelEvent(event)

    # ------------------------------------------------------------ drag drop
    def dragEnterEvent(self, event):
        if event.mimeData().hasFormat(BLOCK_MIME):
            event.acceptProposedAction()
            return

        super().dragEnterEvent(event)

    def dragMoveEvent(self, event):
        if event.mimeData().hasFormat(BLOCK_MIME):
            event.acceptProposedAction()
            return

        super().dragMoveEvent(event)

    def dropEvent(self, event):
        if not event.mimeData().hasFormat(BLOCK_MIME):
            super().dropEvent(event)
            return

        blockKey = bytes(event.mimeData().data(BLOCK_MIME)).decode()
        self.addBlock(blockKey, self.mapToScene(event.position().toPoint()))
        event.acceptProposedAction()
