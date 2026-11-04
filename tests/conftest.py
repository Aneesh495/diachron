"""Make the tests directory importable so `oracle` resolves without packaging."""

import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
