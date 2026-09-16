"""Channel cross-section and hydrodynamic focusing geometry."""

import math


def channel_area(width: float, depth: float) -> float:
    """Cross-section area of a rectangular channel [m^2]."""
    return width * depth


def hydraulic_diameter(width: float, depth: float) -> float:
    """Hydraulic diameter D_h = 4A/P of a rectangular channel [m]."""
    return 4.0 * width * depth / (2.0 * (width + depth))


def core_area(q_sample: float, q_sheath: float, area: float, velocity_ratio: float = 1.0) -> float:
    """Area of the focused sample core [m^2].

    Mass conservation: q_sample = A_core * v_core, with v_core = velocity_ratio * v_mean
    and v_mean = (q_sample + q_sheath) / area.

    velocity_ratio = 1 is plug flow (core occupies the same fraction of the
    area as of the flow); 2 is the centreline of fully developed laminar flow
    in a circular duct, which halves the core area.
    """
    return area * q_sample / (q_sample + q_sheath) / velocity_ratio


def core_diameter(q_sample: float, q_sheath: float, area: float, velocity_ratio: float = 1.0) -> float:
    """Equivalent circular diameter of the focused core [m]."""
    return math.sqrt(4.0 * core_area(q_sample, q_sheath, area, velocity_ratio) / math.pi)
