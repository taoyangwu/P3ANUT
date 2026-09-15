#!/usr/bin/env python3
"""
P3ANUT launcher.

Starts the PyQt node graph interface:

    python P3ANUT.py [pipeline.p3g.yaml]
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))

from p3anut_ui import main

if __name__ == "__main__":
    sys.exit(main())
