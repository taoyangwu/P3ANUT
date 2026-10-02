"""
The graphics items that make up the canvas: blocks, their ports and the
connections between them.

A port is drawn as a regular polygon whose side count comes from its connection
type, so how far the data has travelled through the pipeline is readable at a
glance. Inputs sit on the left of a block and outputs on the right.
"""

import math

from PyQt6.QtCore import QPointF, QRectF, Qt
from PyQt6.QtGui import QBrush, QColor, QFont, QPainterPath, QPen, QPolygonF
from PyQt6.QtWidgets import (
    QGraphicsEllipseItem,
    QGraphicsItem,
    QGraphicsPathItem,
    QGraphicsSimpleTextItem,
)

from . import connections as C

NODE_WIDTH = 250
HEADER_HEIGHT = 28
PORT_SPACING = 24
PORT_RADIUS = 7
BODY_PADDING = 10
SUBTITLE_HEIGHT = 14
BUTTON_STRIP_HEIGHT = 22

_HEADER_COLORS = {
    "fileInput": "#1f77b4",
    "literal": "#ff7f0e",
    "fasta": "#2ca02c",
    "pairedAssembler": "#2ca02c",
    "sequenceCounter": "#9467bd",
    "runUnifier": "#8c564b",
    "volcanoPlot": "#17becf",
    "upsetPlot": "#17becf",
    "rankingPlot": "#17becf",
    "output": "#bcbd22",
}


def portShapePath(sides, radius=PORT_RADIUS):
    """
    A regular polygon with the given number of sides, or a circle when the
    count is zero. Two sides are never used - the shape would read as a line.
    """
    path = QPainterPath()

    if sides < 3:
        path.addEllipse(QPointF(0, 0), radius, radius)
        return path

    #Point the first vertex up so a triangle reads as a triangle
    polygon = QPolygonF()
    for index in range(sides):
        angle = -math.pi / 2 + (2 * math.pi * index / sides)
        polygon.append(QPointF(radius * math.cos(angle), radius * math.sin(angle)))

    path.addPolygon(polygon)
    path.closeSubpath()
    return path


class PortItem(QGraphicsPathItem):
    """One input or output on a block."""

    def __init__(self, nodeItem, spec, isInput, index):
        super().__init__(nodeItem)
        self.nodeItem = nodeItem
        self.spec = spec
        self.isInput = isInput
        self.index = index
        self.edges = []

        self.setAcceptHoverEvents(True)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setZValue(2)

        self.label = QGraphicsSimpleTextItem(spec.label, nodeItem)
        font = QFont()
        font.setPointSize(8)
        self.label.setFont(font)
        self.label.setBrush(QBrush(QColor("#2b2b2b")))

        self.refresh()

    def connectionType(self):
        """
        The type this port actually carries. An Output block's ports are
        declared as ANY and take on whatever is wired into them.
        """
        if self.spec.type is not C.ANY or self.isInput:
            return self.spec.type

        graph = self.nodeItem.canvas.graph
        return graph.resolvedOutputType(self.nodeItem.node.id, self.spec.key)

    def refresh(self):
        connectionType = self.connectionType()

        self.setPath(portShapePath(connectionType.sides))
        self.setBrush(QBrush(QColor(connectionType.color)))
        self.setPen(QPen(QColor("#33000000"), 1))
        self.setToolTip(f"{self.spec.label} - {connectionType.label}")

        y = (self.nodeItem.portsTop() + self.index * PORT_SPACING + PORT_SPACING / 2)

        if self.isInput:
            self.setPos(0, y)
            self.label.setPos(PORT_RADIUS + 6, y - 7)
        else:
            self.setPos(NODE_WIDTH, y)
            width = self.label.boundingRect().width()
            self.label.setPos(NODE_WIDTH - width - PORT_RADIUS - 6, y - 7)

    def scenePos(self):
        return self.mapToScene(QPointF(0, 0))

    def hoverEnterEvent(self, event):
        self.setPen(QPen(QColor("#222222"), 2))
        super().hoverEnterEvent(event)

    def hoverLeaveEvent(self, event):
        self.setPen(QPen(QColor("#33000000"), 1))
        super().hoverLeaveEvent(event)

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.nodeItem.canvas.beginConnection(self)
            event.accept()
            return

        super().mousePressEvent(event)


class ButtonItem(QGraphicsEllipseItem):
    """A small round button drawn on a block, such as the port + and -."""

    def __init__(self, nodeItem, glyph, callback, radius=8):
        super().__init__(-radius, -radius, radius * 2, radius * 2, nodeItem)
        self.callback = callback

        self.setBrush(QBrush(QColor("#ffffff")))
        self.setPen(QPen(QColor("#888888"), 1))
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setZValue(3)

        self.glyph = QGraphicsSimpleTextItem(glyph, self)
        font = QFont()
        font.setPointSize(9)
        font.setBold(True)
        self.glyph.setFont(font)
        rect = self.glyph.boundingRect()
        self.glyph.setPos(-rect.width() / 2, -rect.height() / 2)

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.callback()
            event.accept()
            return

        super().mousePressEvent(event)


class NodeItem(QGraphicsItem):
    """A block on the canvas."""

    def __init__(self, canvas, node):
        super().__init__()
        self.canvas = canvas
        self.node = node
        self.inputPorts = []
        self.outputPorts = []
        self.buttons = []

        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable, True)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable, True)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemSendsGeometryChanges, True)
        self.setPos(node.x, node.y)
        self.setZValue(1)

        self.title = QGraphicsSimpleTextItem(node.label, self)
        titleFont = QFont()
        titleFont.setPointSize(9)
        titleFont.setBold(True)
        self.title.setFont(titleFont)
        self.title.setBrush(QBrush(QColor("#ffffff")))
        self.title.setPos(10, 7)

        self.subtitle = QGraphicsSimpleTextItem("", self)
        subFont = QFont()
        subFont.setPointSize(7)
        self.subtitle.setFont(subFont)
        self.subtitle.setBrush(QBrush(QColor("#555555")))

        self.rebuildPorts()

    # ------------------------------------------------------------- geometry
    def rowCount(self):
        return max(len(self.node.inputs()), len(self.node.outputs()))

    def hasSubtitle(self):
        return bool(self.subtitle.text())

    def portsTop(self):
        """Where the first port row starts, below the header and any subtitle."""
        return (HEADER_HEIGHT + BODY_PADDING
                + (SUBTITLE_HEIGHT if self.hasSubtitle() else 0))

    def height(self):
        #Only the blocks that carry the + and - buttons need room for them
        buttons = BUTTON_STRIP_HEIGHT if self.node.definition.dynamicInputs else 0

        return (self.portsTop() + self.rowCount() * PORT_SPACING
                + BODY_PADDING + buttons)

    def boundingRect(self):
        return QRectF(-PORT_RADIUS - 2, -2,
                      NODE_WIDTH + PORT_RADIUS * 2 + 4, self.height() + 4)

    # ---------------------------------------------------------------- ports
    def rebuildPorts(self):
        for port in self.inputPorts + self.outputPorts:
            port.label.setParentItem(None)
            port.setParentItem(None)
            if port.scene():
                port.scene().removeItem(port)

        for button in self.buttons:
            if button.scene():
                button.scene().removeItem(button)
            button.setParentItem(None)

        #The subtitle shifts the port rows down, so it has to be settled first
        self._refreshSubtitle()

        self.inputPorts = [PortItem(self, spec, True, i)
                           for i, spec in enumerate(self.node.inputs())]
        self.outputPorts = [PortItem(self, spec, False, i)
                            for i, spec in enumerate(self.node.outputs())]
        self.buttons = []

        if self.node.definition.dynamicInputs:
            self._addPortButtons()

        self.prepareGeometryChange()
        self.update()

    def _addPortButtons(self):
        y = self.height() - BUTTON_STRIP_HEIGHT / 2 - 2

        add = ButtonItem(self, "+", lambda: self.canvas.changeInputCount(self, +1))
        add.setPos(NODE_WIDTH - 46, y)
        remove = ButtonItem(self, "−", lambda: self.canvas.changeInputCount(self, -1))
        remove.setPos(NODE_WIDTH - 22, y)

        label = self.node.definition.dynamicLabel.lower()
        add.setToolTip(f"Add another {label} input")
        remove.setToolTip(f"Remove the last {label} input")

        self.buttons = [add, remove]

    def refreshLayout(self):
        """
        Re-apply the layout after a parameter change, keeping the existing port
        items so the edges attached to them stay valid.
        """
        self._refreshSubtitle()

        for port in self.inputPorts + self.outputPorts:
            port.refresh()

        if self.buttons:
            y = self.height() - BUTTON_STRIP_HEIGHT / 2 - 2
            self.buttons[0].setPos(NODE_WIDTH - 46, y)
            self.buttons[1].setPos(NODE_WIDTH - 22, y)

        self.prepareGeometryChange()
        self.update()

    def _refreshSubtitle(self):
        """A one line hint under the title: the chosen file or destination."""
        node = self.node
        text = ""

        if node.blockKey == "fileInput":
            path = node.params.get("filePath")
            text = path.split("/")[-1] if path else "no file selected"
        elif node.blockKey == "output":
            directory = node.params.get("outputDirectory")
            base = node.params.get("baseFileName") or "output"
            text = f"{directory.rstrip('/').split('/')[-1]}/{base}" if directory else "no folder chosen"
        elif node.blockKey == "literal":
            text = f"{node.params.get('valueType', 'String')} = {node.params.get('value', '')}"

        self.subtitle.setText(text)
        self.subtitle.setPos(10, HEADER_HEIGHT + 2)
        self.subtitle.setVisible(bool(text))

    def refreshPortTypes(self):
        for port in self.inputPorts + self.outputPorts:
            port.refresh()

    def allPorts(self):
        return self.inputPorts + self.outputPorts

    # -------------------------------------------------------------- painting
    def paint(self, painter, option, widget=None):
        painter.setRenderHint(painter.RenderHint.Antialiasing, True)

        body = QRectF(0, 0, NODE_WIDTH, self.height())

        path = QPainterPath()
        path.addRoundedRect(body, 8, 8)

        painter.setBrush(QBrush(QColor("#fbfbfb")))
        painter.setPen(QPen(QColor("#3c78d8") if self.isSelected() else QColor("#c8c8c8"),
                            2 if self.isSelected() else 1))
        painter.drawPath(path)

        header = QPainterPath()
        header.addRoundedRect(QRectF(0, 0, NODE_WIDTH, HEADER_HEIGHT + 8), 8, 8)
        painter.setClipRect(QRectF(0, 0, NODE_WIDTH, HEADER_HEIGHT))
        painter.setBrush(QBrush(QColor(_HEADER_COLORS.get(self.node.blockKey, "#666666"))))
        painter.setPen(Qt.PenStyle.NoPen)
        painter.drawPath(header)
        painter.setClipping(False)

    # ---------------------------------------------------------- interaction
    def itemChange(self, change, value):
        if change == QGraphicsItem.GraphicsItemChange.ItemPositionHasChanged:
            self.node.x = self.pos().x()
            self.node.y = self.pos().y()
            self.canvas.refreshEdges()

        return super().itemChange(change, value)

    def mouseDoubleClickEvent(self, event):
        self.canvas.openParameters(self)
        event.accept()

    def mouseReleaseEvent(self, event):
        #A click that did not drag the block opens its parameters
        if (event.button() == Qt.MouseButton.LeftButton
                and (event.screenPos() - event.buttonDownScreenPos(Qt.MouseButton.LeftButton))
                .manhattanLength() < 4):
            self.canvas.openParameters(self)

        super().mouseReleaseEvent(event)


class EdgeItem(QGraphicsPathItem):
    """A connection drawn between two ports."""

    def __init__(self, canvas, edge, sourcePort, targetPort):
        super().__init__()
        self.canvas = canvas
        self.edge = edge
        self.sourcePort = sourcePort
        self.targetPort = targetPort

        self.setZValue(0)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable, True)
        self.setAcceptHoverEvents(True)
        self.setCursor(Qt.CursorShape.PointingHandCursor)

        sourcePort.edges.append(self)
        targetPort.edges.append(self)

        self.refresh()

    def refresh(self):
        start = self.sourcePort.scenePos()
        end = self.targetPort.scenePos()

        #A horizontal bezier, so the curve always leaves right and arrives left
        reach = max(40.0, abs(end.x() - start.x()) * 0.5)

        path = QPainterPath(start)
        path.cubicTo(QPointF(start.x() + reach, start.y()),
                     QPointF(end.x() - reach, end.y()),
                     end)
        self.setPath(path)

        color = QColor(self.sourcePort.connectionType().color)
        width = 3 if self.isSelected() else 2
        self.setPen(QPen(color, width, Qt.PenStyle.SolidLine,
                         Qt.PenCapStyle.RoundCap))
        self.setToolTip(f"{self.sourcePort.connectionType().label} connection - "
                        f"select and press Delete to remove")

    def hoverEnterEvent(self, event):
        pen = self.pen()
        pen.setWidth(4)
        self.setPen(pen)
        super().hoverEnterEvent(event)

    def hoverLeaveEvent(self, event):
        self.refresh()
        super().hoverLeaveEvent(event)
