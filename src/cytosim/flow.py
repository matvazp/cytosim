"""Bulk flow quantities in the channel."""

# Centreline / mean velocity ratio for each supported profile.
# "parabolic" uses the circular-duct value 2.0; a square duct is ~2.1.
VELOCITY_RATIO = {"plug": 1.0, "parabolic": 2.0}


def mean_velocity(q_total: float, area: float) -> float:
    """Mean flow velocity [m/s]."""
    return q_total / area


def reynolds(velocity: float, hydraulic_diameter: float, density: float, viscosity: float) -> float:
    """Reynolds number; laminar (and thus stable focusing) for Re < ~2000."""
    return density * velocity * hydraulic_diameter / viscosity


def velocity_ratio(profile: str) -> float:
    """v_core / v_mean for the given profile name."""
    try:
        return VELOCITY_RATIO[profile]
    except KeyError:
        raise ValueError(f"unknown profile {profile!r}; choose from {list(VELOCITY_RATIO)}") from None


def core_velocity(q_total: float, area: float, profile: str = "plug") -> float:
    """Velocity seen by particles on the channel axis [m/s]."""
    return velocity_ratio(profile) * mean_velocity(q_total, area)
