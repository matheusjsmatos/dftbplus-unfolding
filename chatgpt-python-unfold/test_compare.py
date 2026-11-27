"""Compatibility entry point for the tests now maintained in Testes/."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from Testes.test_spectral_maps import MapChecks

__test__ = False  # pytest collects the canonical module in Testes/.


if __name__ == "__main__":
    unittest.main()
