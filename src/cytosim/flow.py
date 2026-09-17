"""Velocity field in the channel.

Coordinates: the flow runs along z; the cross-section is (x, y) with the
origin on the channel axis, x in [-width/2, width/2] and y in [-depth/2, depth/2].
x is the direction along the laser beam, y the direction across it (the
direction the 100 µm spot width spans).

Two profiles are provided, both normalised so that the mean over the
cross-section equals the bulk velocity q_total / area:

* `RectangularDuct` - fully developed laminar (Poiseuille) flow in a
  rectangular duct, the classic series solution (White, Viscous Fluid Flow,
  eq. 3.48). The fluid sticks to the walls, so the velocity is zero there and
  peaks on the axis at centreline_ratio x the mean (2.096 for a square duct,
  2.0 for a round pipe, 1.5 for parallel plates).
* `PlugFlow` - uniform velocity. Not physical in a laminar duct, but a useful
  lower bound on the particle velocity because the instrument runs at
  Re ~ 2500, where the profile starts to flatten towards turbulent plug flow.
"""

import math

import numpy as np


def mean_velocity(q_total: float, area: float) -> float:
    """Bulk (mean) flow velocity [m/s]."""
    return q_total / area


def reynolds(velocity: float, hydraulic_diameter: float, density: float, viscosity: float) -> float:
    """Reynolds number; laminar (and thus stable focusing) for Re < ~2000."""
    return density * velocity * hydraulic_diameter / viscosity


class PlugFlow:
    """Uniform velocity across the whole cross-section (u / u_mean = 1)."""

    def __init__(self, width: float, depth: float):
        self.width = width
        self.depth = depth
        self.centreline_ratio = 1.0

    def shape(self, x, y) -> np.ndarray:
        """u / u_mean at (x, y) [-]."""
        return np.ones(np.broadcast(np.asarray(x, float), np.asarray(y, float)).shape)

    def velocity(self, x, y, q_total: float) -> np.ndarray:
        """Velocity at (x, y) for a total flow q_total [m/s]."""
        return mean_velocity(q_total, self.width * self.depth) * self.shape(x, y)


class RectangularDuct:
    """Fully developed laminar flow in a rectangular duct (series solution).

    With half-widths a = width/2, b = depth/2 the velocity is

        u(x, y) ∝ Σ_{i odd} (-1)^((i-1)/2) / i^3
                    * [1 - cosh(i π y / 2a) / cosh(i π b / 2a)] * cos(i π x / 2a)

    Each term already vanishes on x = ±a (the cosine) and the bracket makes it
    vanish on y = ±b; the sum converges as 1/i^3 so ~50 terms are ample. The
    mean over the section has a closed form, which is used to normalise:

        <g> = π^3/48 * [1 - 192 a / (π^5 b) * Σ tanh(i π b / 2a) / i^5]
    """

    def __init__(self, width: float, depth: float, n_terms: int = 50):
        self.width = width
        self.depth = depth
        self.n_terms = n_terms
        a, b = width / 2.0, depth / 2.0
        self._i = 2 * np.arange(n_terms) + 1              # 1, 3, 5, ...
        self._k = self._i * math.pi / (2.0 * a)
        self._sign = np.where((self._i // 2) % 2, -1.0, 1.0)   # (-1)^((i-1)/2)
        self._mean = math.pi**3 / 48.0 * (
            1.0 - 192.0 * a / (math.pi**5 * b) * np.sum(np.tanh(self._k * b) / self._i**5)
        )
        self.centreline_ratio = float(self.shape(0.0, 0.0))

    def shape(self, x, y) -> np.ndarray:
        """u / u_mean at (x, y) [-]. Broadcasts over array inputs."""
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        b = self.depth / 2.0
        ay = np.abs(y)
        g = np.zeros(np.broadcast(x, y).shape)
        for i, k, sign in zip(self._i, self._k, self._sign):
            # cosh(k y) / cosh(k b) written with exponentials so large i cannot overflow
            ratio = np.exp(k * (ay - b)) * (1.0 + np.exp(-2.0 * k * ay)) / (1.0 + np.exp(-2.0 * k * b))
            g = g + sign / i**3 * (1.0 - ratio) * np.cos(k * x)
        return g / self._mean

    def velocity(self, x, y, q_total: float) -> np.ndarray:
        """Velocity at (x, y) for a total flow q_total [m/s]."""
        return mean_velocity(q_total, self.width * self.depth) * self.shape(x, y)
