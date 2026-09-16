"""cytosim - flow cytometer channel model for ADC sample-rate estimation.

Pure NumPy/SciPy so it runs both natively and in the browser (Pyodide/stlite).
"""

from .model import Results, simulate
from .params import Params
from . import units

__all__ = ["Params", "Results", "simulate", "units"]
__version__ = "0.1.0"
