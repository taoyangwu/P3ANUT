"""
P3ANUT's PyQt user interface: a single page where the whole pipeline is built
by dragging blocks onto a canvas and wiring them together.
"""

__all__ = ["main"]


def main(argv=None):
    """Launch the application."""
    import sys

    from PyQt6.QtWidgets import QApplication

    from .window import MainWindow

    argv = sys.argv if argv is None else argv

    app = QApplication(argv)
    app.setApplicationName("P3ANUT")

    #An optional pipeline file can be opened straight from the command line
    graphPath = argv[1] if len(argv) > 1 else None

    window = MainWindow(graphPath)
    window.show()

    return app.exec()
