import math

import numpy as np
import pytest

from cytosim import Params, simulate, units
from cytosim import flow, geometry, particle, signal
from cytosim.params import MEASURED_CORE_WIDTH_RANGE


# --- velocity field ----------------------------------------------------------

def _section_mean(profile, n=401):
    xs = np.linspace(-profile.width / 2, profile.width / 2, n)
    ys = np.linspace(-profile.depth / 2, profile.depth / 2, n)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    g = profile.shape(X, Y)
    return np.trapezoid(np.trapezoid(g, ys, axis=1), xs) / (profile.width * profile.depth)


def test_duct_profile_integrates_to_total_flow():
    duct = flow.RectangularDuct(units.um(200), units.um(200))
    assert _section_mean(duct) == pytest.approx(1.0, rel=1e-3)
    q = units.ml_per_min(20)
    area = geometry.channel_area(duct.width, duct.depth)
    assert _section_mean(duct) * flow.mean_velocity(q, area) * area == pytest.approx(q, rel=1e-3)


def test_duct_profile_integrates_for_non_square_section():
    duct = flow.RectangularDuct(units.um(100), units.um(400))
    assert _section_mean(duct) == pytest.approx(1.0, rel=1e-3)


def test_square_duct_centreline_ratio():
    """Literature value for a square duct is 2.096; a round pipe would be 2.0."""
    duct = flow.RectangularDuct(units.um(200), units.um(200))
    assert duct.centreline_ratio == pytest.approx(2.10, abs=0.01)
    assert duct.centreline_ratio == pytest.approx(2.0963, abs=1e-3)


def test_wide_duct_tends_to_parallel_plates():
    """Aspect ratio -> infinity recovers the 1.5 of plane Poiseuille flow."""
    duct = flow.RectangularDuct(units.um(200), units.um(200) * 50)
    assert duct.centreline_ratio == pytest.approx(1.5, rel=0.02)


def test_duct_profile_no_slip_and_symmetry():
    duct = flow.RectangularDuct(units.um(200), units.um(300))
    a, b = duct.width / 2, duct.depth / 2
    assert abs(duct.shape(a, 0.0)) < 1e-9
    assert abs(duct.shape(0.0, b)) < 1e-9
    assert abs(duct.shape(a, b)) < 1e-9
    assert duct.shape(0.3 * a, 0.4 * b) == pytest.approx(duct.shape(-0.3 * a, -0.4 * b))
    assert duct.shape(0.0, 0.0) > duct.shape(0.5 * a, 0.0) > 0.0


def test_plug_flow_is_uniform():
    plug = flow.PlugFlow(units.um(200), units.um(200))
    assert plug.centreline_ratio == 1.0
    assert np.all(plug.shape(np.linspace(-1e-4, 1e-4, 5), 0.0) == 1.0)


# --- focused core --------------------------------------------------------------

def test_core_carries_the_sample_flow():
    duct = flow.RectangularDuct(units.um(200), units.um(200))
    q_s, q_t = units.ul_per_s(10), units.ul_per_s(10) + units.ml_per_min(20)
    rx, ry = geometry.core_size(duct, q_s, q_t)
    assert rx == pytest.approx(ry)
    assert geometry.core_flux(duct, q_t, rx, ry) == pytest.approx(q_s, rel=1e-9)
    rx, ry = geometry.core_size(duct, q_s, q_t, aspect=3.0)
    assert rx / ry == pytest.approx(3.0)
    assert geometry.core_flux(duct, q_t, rx, ry) == pytest.approx(q_s, rel=1e-9)


def test_plug_core_area_equals_flow_fraction():
    """Under plug flow the core occupies the same fraction of the area as of the flow."""
    plug = flow.PlugFlow(units.um(250), units.um(250))
    area = geometry.channel_area(plug.width, plug.depth)
    rx, ry = geometry.core_size(plug, units.ul_per_min(5), units.ul_per_min(500))
    assert math.pi * rx * ry == pytest.approx(area * 0.01, rel=1e-9)


@pytest.mark.parametrize("preset", [Params.eb, Params.pr2])
def test_core_width_matches_instrument(preset):
    """Model core width must land in the observed 10-30 µm for both sample modes.

    No tunable profile flag: the laminar duct field alone gets it right.
    """
    lo, hi = MEASURED_CORE_WIDTH_RANGE
    r = simulate(preset())
    assert lo <= r.core_width <= hi
    assert lo <= r.core_diameter <= hi


def test_plug_flow_overestimates_core_width():
    """Plug flow (the transitional lower bound) misses the measured core width."""
    r = simulate(Params())
    assert r.core_diameter_plug > MEASURED_CORE_WIDTH_RANGE[1]
    assert r.core_diameter_plug > r.core_diameter


def test_core_velocity_spread_and_bounds():
    r = simulate(Params())
    assert r.velocity == pytest.approx(r.velocity_ratio * r.velocity_bulk)
    assert r.velocity_bulk < r.velocity_min < r.velocity_mean < r.velocity
    # a core this small sits on the flat top of the profile: spread of a few %
    assert (r.velocity - r.velocity_min) / r.velocity < 0.05
    assert r.fwhm < r.fwhm_max
    assert r.f_s_required_plug < r.f_s_required


def test_ribbon_core_is_thin_across_the_beam():
    r = simulate(Params(core_aspect=4.0))
    assert r.core_height == pytest.approx(4.0 * r.core_width)
    assert r.core_diameter == pytest.approx(simulate(Params()).core_diameter, rel=0.05)


def test_impossible_core_raises():
    duct = flow.RectangularDuct(units.um(200), units.um(200))
    with pytest.raises(ValueError):
        geometry.core_size(duct, units.ml_per_min(20), units.ml_per_min(20))


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
    assert r.velocity_ratio == pytest.approx(2.10, abs=0.01)
    # The instrument runs at Re ~2500: transitional, not fully turbulent (>4000).
    assert 2000 < r.reynolds < 4000
    assert 1e6 < r.f_s_required < 1e8   # MHz..tens of MHz
    assert r.pulse.max() == pytest.approx(1.0)
    assert set(r.summary()) == set(r.__dataclass_fields__) - {"t", "pulse"}
