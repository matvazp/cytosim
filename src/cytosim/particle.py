"""Particle transit through the laser spot and the resulting detector pulse.

A particle of diameter d moving at velocity v along y through a Gaussian
beam of 1/e^2 radius w produces a signal proportional to the beam intensity
integrated over the particle extent (top-hat particle, uniform emission):

    S(t) ∝ ∫_{vt-d/2}^{vt+d/2} exp(-2 y^2 / w^2) dy
         = w sqrt(pi/8) [erf(sqrt2 (vt + d/2)/w) - erf(sqrt2 (vt - d/2)/w)]

i.e. a top-hat convolved with a Gaussian. For d << w the pulse is Gaussian
with sigma_t = w / (2 v); for d >> w it is a top-hat of duration d / v.
"""

import math

import numpy as np
from scipy.special import erf

from .optics import waist_radius


def gaussian_sigma_t(spot_diameter: float, velocity: float) -> float:
    """Pulse sigma [s] in the point-particle limit."""
    return waist_radius(spot_diameter) / (2.0 * velocity)


def transit_time(spot_diameter: float, particle_diameter: float, velocity: float) -> float:
    """Simple 'edge-to-edge' transit time [s]: (spot + particle) / v."""
    return (spot_diameter + particle_diameter) / velocity


def pulse_shape(t: np.ndarray, velocity: float, particle_diameter: float, spot_diameter: float) -> np.ndarray:
    """Normalised detector pulse (peak = 1) as a function of time [s].

    t = 0 corresponds to the particle centre on the beam axis.
    """
    w = waist_radius(spot_diameter)
    d = particle_diameter
    k = math.sqrt(2.0) / w
    y = velocity * np.asarray(t, dtype=float)
    s = erf(k * (y + d / 2.0)) - erf(k * (y - d / 2.0))
    peak = 2.0 * erf(k * d / 2.0)
    return s / peak


def pulse_time_axis(velocity: float, particle_diameter: float, spot_diameter: float,
                    n: int = 2001, span: float = 3.0) -> np.ndarray:
    """Symmetric time grid covering ±span x the transit time."""
    half = span * transit_time(spot_diameter, particle_diameter, velocity) / 2.0
    return np.linspace(-half, half, n)


def pulse_fwhm(t: np.ndarray, s: np.ndarray) -> float:
    """Full width at half maximum of a sampled unimodal pulse [s]."""
    half = 0.5 * s.max()
    above = np.flatnonzero(s >= half)
    i0, i1 = above[0], above[-1]

    def _cross(i, j):  # linear interpolation of the half-max crossing between samples i and j
        return t[i] + (half - s[i]) * (t[j] - t[i]) / (s[j] - s[i])

    t_left = _cross(i0 - 1, i0) if i0 > 0 else t[i0]
    t_right = _cross(i1, i1 + 1) if i1 < len(s) - 1 else t[i1]
    return t_right - t_left
