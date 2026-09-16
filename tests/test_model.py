import math

import numpy as np
import pytest

from cytosim import Params, simulate, units
from cytosim import geometry, particle, signal
from cytosim.params import MEASURED_CORE_WIDTH_RANGE


def test_core_area_conserves_flow_fraction():
    area = geometry.channel_area(units.um(250), units.um(250))
    core = geometry.core_area(units.ul_per_min(5), units.ul_per_min(495), area)
    assert core == pytest.approx(area * 0.01)
    # a core moving twice as fast carries the same flow through half the area
    assert geometry.core_area(units.ul_per_min(5), units.ul_per_min(495), area, 2.0) == pytest.approx(area * 0.005)


@pytest.mark.parametrize("preset", [Params.eb, Params.pr2])
def test_core_width_matches_instrument(preset):
    """Model core width must land in the observed 10-30 µm for both sample modes."""
    lo, hi = MEASURED_CORE_WIDTH_RANGE
    r = simulate(preset())
    assert lo <= r.core_diameter <= hi


def test_plug_flow_overestimates_core_width():
    """Documents why the default profile is parabolic: plug flow misses the measurement."""
    r = simulate(Params(velocity_profile="plug"))
    assert r.core_diameter > MEASURED_CORE_WIDTH_RANGE[1]


def test_unknown_profile_raises():
    with pytest.raises(ValueError):
        simulate(Params(velocity_profile="turbulent"))


def test_point_particle_pulse_is_gaussian():
    v, spot = 0.1, units.um(10)
    t = np.linspace(-1e-4, 1e-4, 4001)
    s = particle.pulse_shape(t, v, units.um(0.01), spot)
    sigma = particle.gaussian_sigma_t(spot, v)
    expected = np.exp(-t**2 / (2 * sigma**2))
    assert np.allclose(s, expected, atol=1e-4)


def test_large_particle_pulse_is_tophat():
    v, spot, d = 0.1, units.um(1), units.um(100)
    t = particle.pulse_time_axis(v, d, spot)
    s = particle.pulse_shape(t, v, d, spot)
    assert particle.pulse_fwhm(t, s) == pytest.approx(d / v, rel=1e-2)


def test_gaussian_bandwidth():
    sigma = 1e-5
    f = signal.gaussian_bandwidth_3db(sigma)
    # power spectrum at f_3dB is exactly one half
    assert math.exp(-4 * math.pi**2 * sigma**2 * f**2) == pytest.approx(0.5)


def test_faster_flow_needs_higher_sample_rate():
    slow = simulate(Params(q_sheath=units.ul_per_min(200)))
    fast = simulate(Params(q_sheath=units.ul_per_min(800)))
    assert fast.f_s_required > slow.f_s_required


def test_defaults_are_plausible():
    r = simulate(Params())
    # The instrument runs at Re ~2500: transitional, not fully turbulent (>4000).
    assert 2000 < r.reynolds < 4000
    assert 1e6 < r.f_s_required < 1e8   # MHz..tens of MHz
    assert r.pulse.max() == pytest.approx(1.0)
    assert set(r.summary()) == set(r.__dataclass_fields__) - {"t", "pulse"}
