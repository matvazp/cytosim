"""Channel cross-section and hydrodynamic focusing geometry.

The sample is injected on the channel axis and the sheath squeezes it into a
thin core. Mass conservation fixes the core size: the flow *through the core
region* must equal the sample flow rate,

    q_sample = ∫∫_core u(x, y) dA.

Because the core sits where the fluid is fastest, the same µL/s needs less
area than the plug-flow estimate A * q_sample / q_total. The core cross-section
is modelled as an ellipse centred on the axis with semi-axes (rx, ry);
aspect = rx / ry = 1 is a circular core, > 1 a ribbon that is tall along the
laser beam (x) and thin across it (y).
"""

import math

import numpy as np
from scipy.optimize import brentq


def channel_area(width: float, depth: float) -> float:
    """Cross-section area of a rectangular channel [m^2]."""
    return width * depth


def hydraulic_diameter(width: float, depth: float) -> float:
    """Hydraulic diameter D_h = 4A/P of a rectangular channel [m]."""
    return 4.0 * width * depth / (2.0 * (width + depth))


def _ellipse_grid(rx: float, ry: float, n_r: int, n_theta: int):
    """Midpoint-rule quadrature over an ellipse: points (x, y) and area weights."""
    rho = (np.arange(n_r) + 0.5) / n_r
    theta = (np.arange(n_theta) + 0.5) / n_theta * 2.0 * math.pi
    R, T = np.meshgrid(rho, theta, indexing="ij")
    x = rx * R * np.cos(T)
    y = ry * R * np.sin(T)
    w = rx * ry * R * (1.0 / n_r) * (2.0 * math.pi / n_theta)
    return x, y, w


def core_flux(profile, q_total: float, rx: float, ry: float, n_r: int = 64, n_theta: int = 32) -> float:
    """Volume flow [m^3/s] through the ellipse with semi-axes (rx, ry)."""
    x, y, w = _ellipse_grid(rx, ry, n_r, n_theta)
    return float(np.sum(profile.velocity(x, y, q_total) * w))


def core_size(profile, q_sample: float, q_total: float, aspect: float = 1.0) -> tuple[float, float]:
    """Core semi-axes (rx, ry) [m] such that the flow through the core is q_sample.

    aspect = rx / ry. Solved by bracketing on rx between 0 and the largest
    ellipse of that aspect ratio that fits in the channel.
    """
    rx_max = min(profile.width / 2.0, aspect * profile.depth / 2.0)

    def residual(rx):
        return core_flux(profile, q_total, rx, rx / aspect) - q_sample

    if residual(rx_max) < 0.0:
        raise ValueError("sample flow exceeds what the channel can carry through a core of this aspect ratio")
    rx = brentq(residual, 1e-12 * rx_max, rx_max, xtol=1e-12)
    return rx, rx / aspect


def core_diameter(rx: float, ry: float) -> float:
    """Equivalent circular diameter of an elliptical core [m]."""
    return 2.0 * math.sqrt(rx * ry)


def core_velocities(profile, q_total: float, rx: float, ry: float,
                    n_r: int = 64, n_theta: int = 32) -> tuple[np.ndarray, np.ndarray]:
    """Velocities of particles crossing the core, with their flux weights.

    Particles are carried with the fluid, so the number of particles passing
    through an area element dA per second is proportional to u dA. The
    returned weights therefore give the velocity distribution of *detected*
    particles; they sum to the core flow rate.
    """
    x, y, w = _ellipse_grid(rx, ry, n_r, n_theta)
    v = profile.velocity(x, y, q_total)
    return v.ravel(), (v * w).ravel()


def core_velocity_extremes(profile, q_total: float, rx: float, ry: float, n_theta: int = 256) -> tuple[float, float]:
    """(v_min, v_max) [m/s] over the core: the boundary minimum and the axis value."""
    theta = np.linspace(0.0, 2.0 * math.pi, n_theta, endpoint=False)
    v_edge = profile.velocity(rx * np.cos(theta), ry * np.sin(theta), q_total)
    v_axis = profile.velocity(0.0, 0.0, q_total)
    return float(np.min(v_edge)), float(v_axis)
