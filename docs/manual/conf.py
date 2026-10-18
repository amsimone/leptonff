from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bindings/python"))
project = "leptonff"
author = "Simone Amoroso"
release = "0.26.0"
extensions = ["myst_parser", "sphinx.ext.mathjax", "sphinx.ext.autodoc"]
myst_enable_extensions = ["dollarmath", "amsmath"]
html_theme = "alabaster"
exclude_patterns = []
nitpicky = False
