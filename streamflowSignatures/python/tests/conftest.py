"""
Pytest configuration and fixtures for streamflow signatures tests.
"""

import sys
from pathlib import Path

# Ensure the package is importable
sys.path.insert(0, str(Path(__file__).parent.parent))
