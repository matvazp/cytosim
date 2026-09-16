"""Unit helpers.

The package works in SI internally (m, s, m^3/s, kg/m^3, Pa.s).
These helpers convert the "human" units used in the notebook and the app
(µm, µL/min, µL/s, mL/min) at the edges.
"""

UM = 1e-6                  # µm -> m
MM = 1e-3                  # mm -> m
UL_PER_MIN = 1e-9 / 60.0   # µL/min -> m^3/s
UL_PER_S = 1e-9            # µL/s   -> m^3/s
ML_PER_MIN = 1e-6 / 60.0   # mL/min -> m^3/s


def um(value: float) -> float:
    """Micrometres to metres."""
    return value * UM


def ul_per_min(value: float) -> float:
    """Microlitres per minute to m^3/s."""
    return value * UL_PER_MIN


def ul_per_s(value: float) -> float:
    """Microlitres per second to m^3/s."""
    return value * UL_PER_S


def ml_per_min(value: float) -> float:
    """Millilitres per minute to m^3/s."""
    return value * ML_PER_MIN


def to_um(value_m: float) -> float:
    return value_m / UM


def to_ul_per_min(value_m3_s: float) -> float:
    return value_m3_s / UL_PER_MIN


def to_ml_per_min(value_m3_s: float) -> float:
    return value_m3_s / ML_PER_MIN
