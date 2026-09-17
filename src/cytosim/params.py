"""Simulation parameters.

One dataclass holds every user-adjustable input, with defaults taken from the
current real-life instrument. The notebook and the stlite app both build their
sliders from this class, so adding a parameter here is the single place to
extend the model.

UNITS: every field is STORED in SI (m, s, m^3/s, kg/m^3, Pa.s, W).
The number you TYPE goes through a converter so it can be in human units:
    um(200.0)          -> 200 µm, stored as 2.0e-4 m
    ul_per_min(600.0)  -> 600 µL/min, stored as 1.0e-8 m^3/s
Never pass a bare number for a length or flow rate: Params(channel_width=200)
would mean a 200 m wide channel and silently produce nonsense.
"""

from dataclasses import dataclass, asdict, replace

from .units import um, ul_per_s, ml_per_min

# --- Reference instrument facts that are not model inputs -------------------

# Observed width of the focused sample stream in the 200x200 µm channel.
# Used to validate the hydrodynamic-focusing model, not fed into it.
MEASURED_CORE_WIDTH_RANGE = (um(10.0), um(30.0))   # m

# Particle populations in raw milk. Diameters in m, concentration per m^3 of
# undiluted milk. Not used by the pulse-timing model yet; kept here for the
# event-rate / coincidence and photon-budget work to come.
MILK_PARTICLES = {
    #                        diameter (typ.)   diameter range      conc. [1/m^3]  (per mL)
    "casein micelle":   dict(d=um(0.2),   d_range=(um(0.1),  um(0.3)),   n=None),
    "fat globule":      dict(d=um(3.0),   d_range=(um(0.1),  um(20.0)),  n=1e9 * 1e6),   # ~1e9/mL
    "somatic cell":     dict(d=um(7.5),   d_range=(um(5.0),  um(10.0)),  n=2e5 * 1e6),   # 5e4 - 2e6/mL
    "bacterium":        dict(d=um(1.0),   d_range=(um(0.5),  um(5.0)),   n=1e5 * 1e6),   # 5e3 - 1e7/mL
}


@dataclass
class Params:
    # --- Channel geometry (rectangular cross-section) ---
    channel_width: float = um(200.0)       # stored m; type µm via um()
    channel_depth: float = um(200.0)       # stored m; type µm via um()

    # --- Flow rates ---
    q_sample: float = ul_per_s(10.0)       # stored m^3/s; diluted sample stream (10 µL/s = 600 µL/min)
    q_sheath: float = ml_per_min(20.0)     # stored m^3/s; sheath at ~1.5 bar
    sample_dilution: float = 2.0           # (milk + buffer) / milk volume; 1:1 EB -> 2.0, 1:3.2 PR2 -> 4.2

    # --- Fluid (sheath-dominated; water at ~37 °C when it reaches the channel) ---
    temperature: float = 37.0              # °C, informative
    density: float = 993.0                 # kg/m^3
    viscosity: float = 0.69e-3             # Pa.s

    # Focused core cross-section: ellipse with aspect = extent along the laser
    # beam (x) / extent across it (y). 1.0 = circular core; > 1 = ribbon that
    # is thin across the beam and tall along it. Injector geometry unknown, so
    # circular is assumed until measured.
    core_aspect: float = 1.0

    # --- Laser spot (1/e^2 intensity diameters of the oval spot) ---
    spot_height: float = um(20.0)          # stored m; along the flow direction
    spot_width: float = um(100.0)          # stored m; across the flow direction
    wavelength: float = 450e-9             # m (not used by the timing model yet)
    laser_power: float = 40e-3             # W, ">40 mW" (not used by the timing model yet)

    # --- Particle under study ---
    particle_diameter: float = um(1.0)     # stored m; see MILK_PARTICLES for typical values

    # --- Sampling criterion ---
    samples_per_fwhm: float = 10.0         # samples across the pulse FWHM
    oversampling: float = 5.0              # multiple of the Nyquist rate for the -3 dB bandwidth

    def as_dict(self) -> dict:
        return asdict(self)

    @property
    def q_total(self) -> float:
        return self.q_sample + self.q_sheath

    # --- Instrument presets ---------------------------------------------------

    @classmethod
    def eb(cls, **overrides) -> "Params":
        """Sample diluted 1:1 with EB, 10 µL/s at 40 °C (the default)."""
        return replace(cls(), q_sample=ul_per_s(10.0), sample_dilution=2.0, **overrides)

    @classmethod
    def pr2(cls, **overrides) -> "Params":
        """Sample diluted 1:3.2 with PR2, 7 µL/s at 40 °C."""
        return replace(cls(), q_sample=ul_per_s(7.0), sample_dilution=4.2, **overrides)
