"""Laser spot description.

The spot is modelled as a Gaussian intensity profile along the flow axis y:
    I(y) = I0 * exp(-2 y^2 / w^2)
where w is the 1/e^2 radius (spot_height = 2w).
"""

import math

import numpy as np


def waist_radius(spot_diameter: float) -> float:
    """1/e^2 radius from the 1/e^2 diameter [m]."""
    return spot_diameter / 2.0


def gaussian_profile(y: np.ndarray, spot_diameter: float) -> np.ndarray:
    """Normalised intensity along the flow axis (peak = 1)."""
    w = waist_radius(spot_diameter)
    return np.exp(-2.0 * y**2 / w**2)


def beam_fwhm(spot_diameter: float) -> float:
    """Full width at half maximum of the intensity profile [m]."""
    w = waist_radius(spot_diameter)
    return math.sqrt(2.0 * math.log(2.0)) * w
