"""
The parameter popup shown when a block is clicked.

Only the parameters the block's own script reads are listed - not the whole
configuration file. The values belong to that one block: editing them here never
affects another block of the same type, nor the YAML template they were seeded
from.
"""

import os

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from .config import coerce


class _BooleanListEdit(QLineEdit):
    """Comma separated true/false values, one per connected file."""

    def __init__(self, value, parent=None):
        super().__init__(parent)
        self.setText(", ".join("true" if v else "false" for v in (value or [])))
        self.setPlaceholderText("empty includes every file - e.g. true, false, true")


class ParameterDialog(QDialog):
    """Edits one block instance's parameters."""

    def __init__(self, node, template, parent=None):
        super().__init__(parent)
        self.node = node
        self.template = template
        self.editors = {}

        self.setWindowTitle(f"{node.label} parameters")
        self.setMinimumWidth(520)

        layout = QVBoxLayout(self)

        if node.definition.description:
            blurb = QLabel(node.definition.description)
            blurb.setWordWrap(True)
            blurb.setStyleSheet("color: #555; padding-bottom: 6px;")
            layout.addWidget(blurb)

        #The file and folder pickers that belong on the block itself
        if node.blockKey == "fileInput":
            layout.addLayout(self._filePicker())
        elif node.blockKey == "output":
            layout.addLayout(self._folderPicker())

        described = template.describe(node.definition.configSection,
                                      node.definition.configKeys)

        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)

        for key, record in described.items():
            #Handled by the dedicated pickers above
            if node.blockKey == "output" and key in ("outputDirectory", "baseFileName"):
                continue

            editor = self._editorFor(key, record)
            self.editors[key] = (editor, record)

            label = QLabel(key)
            label.setToolTip(record.get("Description", ""))
            editor.setToolTip(record.get("Description", ""))
            form.addRow(label, editor)

        if not described:
            form.addRow(QLabel("This block has no configurable parameters."))

        container = QWidget()
        container.setLayout(form)

        scroll = QScrollArea()
        scroll.setWidget(container)
        scroll.setWidgetResizable(True)
        layout.addWidget(scroll, 1)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok
            | QDialogButtonBox.StandardButton.Cancel
            | QDialogButtonBox.StandardButton.RestoreDefaults)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        buttons.button(QDialogButtonBox.StandardButton.RestoreDefaults).clicked.connect(
            self._restoreDefaults)
        layout.addWidget(buttons)

    # ----------------------------------------------------------- the pickers
    def _filePicker(self):
        row = QHBoxLayout()
        self.fileEdit = QLineEdit(self.node.params.get("filePath", ""))
        self.fileEdit.setPlaceholderText("no file selected")

        browse = QPushButton("Choose file")
        browse.clicked.connect(self._chooseFile)

        row.addWidget(QLabel("File"))
        row.addWidget(self.fileEdit, 1)
        row.addWidget(browse)
        return row

    def _chooseFile(self):
        start = self.node.params.get("startingDirectory") or os.path.expanduser("~")
        path, _ = QFileDialog.getOpenFileName(self, "Select an input file", start)
        if path:
            self.fileEdit.setText(path)

    def _folderPicker(self):
        outer = QVBoxLayout()

        row = QHBoxLayout()
        self.folderEdit = QLineEdit(self.node.params.get("outputDirectory", ""))
        self.folderEdit.setPlaceholderText("no folder chosen")

        browse = QPushButton("\U0001F4C1 Choose folder")
        browse.clicked.connect(self._chooseFolder)

        row.addWidget(QLabel("Folder"))
        row.addWidget(self.folderEdit, 1)
        row.addWidget(browse)
        outer.addLayout(row)

        nameRow = QHBoxLayout()
        self.baseNameEdit = QLineEdit(self.node.params.get("baseFileName", "output"))
        nameRow.addWidget(QLabel("Base file name"))
        nameRow.addWidget(self.baseNameEdit, 1)
        outer.addLayout(nameRow)

        note = QLabel("The extension is set automatically from the connection "
                      "feeding this block.")
        note.setStyleSheet("color: #777; font-size: 11px;")
        outer.addWidget(note)

        return outer

    def _chooseFolder(self):
        start = self.folderEdit.text() or os.path.expanduser("~")
        path = QFileDialog.getExistingDirectory(self, "Select a destination folder", start)
        if path:
            self.folderEdit.setText(path)

    # ---------------------------------------------------------- the editors
    def _editorFor(self, key, record):
        declared = record.get("Type", "String")
        current = self.node.params.get(key, record.get("Value"))

        if declared == "Boolean":
            editor = QCheckBox()
            editor.setChecked(bool(current))
            return editor

        if declared == "Choice":
            editor = QComboBox()
            options = [str(o) for o in record.get("Options", [])]
            editor.addItems(options)
            if str(current) in options:
                editor.setCurrentText(str(current))
            return editor

        if declared == "Int":
            editor = QSpinBox()
            editor.setRange(-1_000_000, 1_000_000)
            editor.setValue(int(current if current is not None else 0))
            return editor

        if declared == "Float":
            editor = QDoubleSpinBox()
            editor.setRange(-1_000_000.0, 1_000_000.0)
            editor.setDecimals(4)
            editor.setValue(float(current if current is not None else 0.0))
            return editor

        if declared == "BooleanList":
            return _BooleanListEdit(current)

        editor = QLineEdit("" if current is None else str(current))
        return editor

    def _restoreDefaults(self):
        """Reset this block back to the template values, leaving others alone."""
        defaults = self.template.defaults(self.node.definition.configSection,
                                          self.node.definition.configKeys)

        for key, (editor, record) in self.editors.items():
            value = defaults.get(key)

            if isinstance(editor, QCheckBox):
                editor.setChecked(bool(value))
            elif isinstance(editor, QComboBox):
                editor.setCurrentText(str(value))
            elif isinstance(editor, QSpinBox):
                editor.setValue(int(value or 0))
            elif isinstance(editor, QDoubleSpinBox):
                editor.setValue(float(value or 0.0))
            elif isinstance(editor, _BooleanListEdit):
                editor.setText(", ".join("true" if v else "false" for v in (value or [])))
            else:
                editor.setText("" if value is None else str(value))

    # ------------------------------------------------------------- applying
    def applyTo(self, node):
        if node.blockKey == "fileInput":
            node.params["filePath"] = self.fileEdit.text().strip()
        elif node.blockKey == "output":
            node.params["outputDirectory"] = self.folderEdit.text().strip()
            node.params["baseFileName"] = self.baseNameEdit.text().strip() or "output"

        for key, (editor, record) in self.editors.items():
            declared = record.get("Type", "String")

            if isinstance(editor, QCheckBox):
                node.params[key] = editor.isChecked()
            elif isinstance(editor, QComboBox):
                node.params[key] = coerce(editor.currentText(), declared)
            elif isinstance(editor, (QSpinBox, QDoubleSpinBox)):
                node.params[key] = editor.value()
            else:
                node.params[key] = coerce(editor.text(), declared)
