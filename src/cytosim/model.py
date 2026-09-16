"""Top-level model: Params in, Results out.

This is the single entry point used by the notebook and the app.
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
    core_diameter: float       # m
    velocity: float            # m/s
    reynolds: float
    # pulse
    transit_time: float        # s
    sigma_t: float             # s
    fwhm: float                # s
    f_3db: float               # Hz
    # sampling
    f_s_fwhm: float            # Hz, from samples-per-FWHM criterion
    f_s_bandwidth: float       # Hz, from oversampled Nyquist criterion
    f_s_required: float        # Hz, max of the two
    # waveform for plotting
    t: np.ndarray              # s
    pulse: np.ndarray          # normalised

    def summary(self) -> dict:
        """Scalars only, for tables and the app."""
        d = asdict(self)
        d.pop("t")
        d.pop("pulse")
        return d


def simulate(p: Params) -> Results:
    area = geometry.channel_area(p.channel_width, p.channel_depth)
    d_h = geometry.hydraulic_diameter(p.channel_width, p.channel_depth)
    ratio = flow.velocity_ratio(p.velocity_profile)
    d_core = geometry.core_diameter(p.q_sample, p.q_sheath, area, ratio)

    v = flow.core_velocity(p.q_total, area, p.velocity_profile)
    re = flow.reynolds(flow.mean_velocity(p.q_total, area), d_h, p.density, p.viscosity)

    t = particle.pulse_time_axis(v, p.particle_diameter, p.spot_height)
    s = particle.pulse_shape(t, v, p.particle_diameter, p.spot_height)
    fwhm = particle.pulse_fwhm(t, s)
    sigma_t = signal.pulse_sigma_t(t, s)
    f_3db = signal.gaussian_bandwidth_3db(sigma_t)

    f_s_fwhm = signal.sample_rate_from_fwhm(fwhm, p.samples_per_fwhm)
    f_s_bw = signal.sample_rate_from_bandwidth(f_3db, p.oversampling)

    return Results(
        area=area,
        hydraulic_diameter=d_h,
        core_diameter=d_core,
        velocity=v,
        reynolds=re,
        transit_time=particle.transit_time(p.spot_height, p.particle_diameter, v),
        sigma_t=sigma_t,
        fwhm=fwhm,
        f_3db=f_3db,
        f_s_fwhm=f_s_fwhm,
        f_s_bandwidth=f_s_bw,
        f_s_required=max(f_s_fwhm, f_s_bw),
        t=t,
        pulse=s,
    )
