"""Top-level model: Params in, Results out.

This is the single entry point used by the notebook and the app.

The velocity field is the laminar rectangular-duct solution (flow.RectangularDuct).
The sample core is sized by mass conservation with that field (geometry.core_size),
and particles inside the core see a range of velocities: fastest on the axis,
slowest at the core edge. The ADC requirement follows from the fastest
particle (shortest pulse). Because the instrument runs at Re ~ 2500
(transitional), the same numbers for plug flow are reported alongside as a
lower bound.
"""

from dataclasses import dataclass, asdict

import numpy as np

from . import flow, geometry, particle, signal
from .params import Params


@dataclass
class Results:
    # geometry / flow
    area: float                # m^2
    hydraulic_diameter: float  # m
    reynolds: float
    velocity_bulk: float       # m/s, q_total / area
    velocity_ratio: float      # centreline / bulk for the duct profile (~2.10 square)
    # focused core (laminar profile)
    core_width: float          # m, extent across the laser beam (y)
    core_height: float         # m, extent along the laser beam (x)
    core_diameter: float       # m, equivalent circular diameter
    core_diameter_plug: float  # m, plug-flow lower-bound estimate
    # particle velocities across the core
    velocity: float            # m/s, fastest particle (channel axis)
    velocity_min: float        # m/s, slowest particle (core edge)
    velocity_mean: float       # m/s, flux-weighted mean over the core
    # pulse of the fastest particle
    transit_time: float        # s
    sigma_t: float             # s
    fwhm: float                # s, shortest pulse (fastest particle)
    fwhm_max: float            # s, longest pulse (slowest particle in the core)
    f_3db: float               # Hz
    # sampling
    f_s_fwhm: float            # Hz, from samples-per-FWHM criterion
    f_s_bandwidth: float       # Hz, from oversampled Nyquist criterion
    f_s_required: float        # Hz, max of the two, fastest particle
    f_s_required_plug: float   # Hz, same criterion under plug flow (lower bound)
    # waveform of the fastest particle, for plotting
    t: np.ndarray              # s
    pulse: np.ndarray          # normalised

    def summary(self) -> dict:
        """Scalars only, for tables and the app."""
        d = asdict(self)
        d.pop("t")
        d.pop("pulse")
        return d


def _pulse(p: Params, v: float):
    """Time axis, normalised pulse and FWHM for a particle moving at v."""
    t = particle.pulse_time_axis(v, p.particle_diameter, p.spot_height)
    s = particle.pulse_shape(t, v, p.particle_diameter, p.spot_height)
    return t, s, particle.pulse_fwhm(t, s)


def _sample_rate(p: Params, fwhm: float, sigma_t: float) -> tuple[float, float, float]:
    f_s_fwhm = signal.sample_rate_from_fwhm(fwhm, p.samples_per_fwhm)
    f_s_bw = signal.sample_rate_from_bandwidth(signal.gaussian_bandwidth_3db(sigma_t), p.oversampling)
    return f_s_fwhm, f_s_bw, max(f_s_fwhm, f_s_bw)


def simulate(p: Params) -> Results:
    area = geometry.channel_area(p.channel_width, p.channel_depth)
    d_h = geometry.hydraulic_diameter(p.channel_width, p.channel_depth)
    v_bulk = flow.mean_velocity(p.q_total, area)
    re = flow.reynolds(v_bulk, d_h, p.density, p.viscosity)

    duct = flow.RectangularDuct(p.channel_width, p.channel_depth)
    rx, ry = geometry.core_size(duct, p.q_sample, p.q_total, p.core_aspect)
    v_min, v_max = geometry.core_velocity_extremes(duct, p.q_total, rx, ry)
    v_core, w_core = geometry.core_velocities(duct, p.q_total, rx, ry)
    v_mean = float(np.sum(v_core * w_core) / np.sum(w_core))

    t, s, fwhm = _pulse(p, v_max)
    sigma_t = signal.pulse_sigma_t(t, s)
    f_3db = signal.gaussian_bandwidth_3db(sigma_t)
    f_s_fwhm, f_s_bw, f_s_req = _sample_rate(p, fwhm, sigma_t)
    _, _, fwhm_max = _pulse(p, v_min)

    # plug-flow lower bound: uniform velocity v_bulk, core from the flow fraction
    plug = flow.PlugFlow(p.channel_width, p.channel_depth)
    rx_plug, ry_plug = geometry.core_size(plug, p.q_sample, p.q_total, p.core_aspect)
    t_plug, s_plug, fwhm_plug = _pulse(p, v_bulk)
    _, _, f_s_plug = _sample_rate(p, fwhm_plug, signal.pulse_sigma_t(t_plug, s_plug))

    return Results(
        area=area,
        hydraulic_diameter=d_h,
        reynolds=re,
        velocity_bulk=v_bulk,
        velocity_ratio=duct.centreline_ratio,
        core_width=2.0 * ry,
        core_height=2.0 * rx,
        core_diameter=geometry.core_diameter(rx, ry),
        core_diameter_plug=geometry.core_diameter(rx_plug, ry_plug),
        velocity=v_max,
        velocity_min=v_min,
        velocity_mean=v_mean,
        transit_time=particle.transit_time(p.spot_height, p.particle_diameter, v_max),
        sigma_t=sigma_t,
        fwhm=fwhm,
        fwhm_max=fwhm_max,
        f_3db=f_3db,
        f_s_fwhm=f_s_fwhm,
        f_s_bandwidth=f_s_bw,
        f_s_required=f_s_req,
        f_s_required_plug=f_s_plug,
        t=t,
        pulse=s,
    )
