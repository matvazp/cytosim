# cytosim — ADC sample-rate model for a milk flow cytometer channel

[![CI](https://github.com/matvazp/cytosim/actions/workflows/ci.yml/badge.svg)](https://github.com/matvazp/cytosim/actions/workflows/ci.yml)

Estimates the ADC sampling rate needed to digitise the detector pulses of a
hydrodynamically focused flow cytometer channel, as a function of channel
geometry, flow rates, laser spot size and particle size.

## Status

**v0.1 — first-pass physics, numbers not yet reviewed.** The package runs and
is tested, but the sampling criterion (10 samples per pulse FWHM, 5× Nyquist
oversampling) is a placeholder to be agreed with signal processing, and the
model does not yet include detector optics, noise or event rate. Current
headline with instrument defaults: ~0.7 µs pulse FWHM for a 1 µm particle,
**~14 MHz** required sample rate. See *Model assumptions* below before quoting.

Roadmap:
1. Agree the sampling criterion; add event rate / coincidence from milk particle concentrations.
2. Detector collection geometry and photon budget (laser power, wavelength are already in `Params`).
3. Browser app (stlite) with an animated channel view, shareable by link.

## Layout

```
src/cytosim/      physics package (pure NumPy/SciPy — runs in Pyodide/stlite)
  params.py       Params dataclass: every user-adjustable input, SI units
  geometry.py     cross-section, hydraulic diameter, focused core size
  flow.py         velocity, Reynolds number
  optics.py       Gaussian laser spot
  particle.py     transit, detector pulse shape (top-hat ⊗ Gaussian)
  signal.py       pulse width/bandwidth → sample-rate criteria
  model.py        simulate(Params) -> Results
tests/            pytest
notebooks/        exploration / derivation notebooks
app/              stlite (Streamlit-in-browser) app — to come
```

## Setup (Windows, PowerShell)

Clone to a short path outside OneDrive (e.g. `C:\Users\<you>\dev\cytosim`).
Long paths hit the Windows 260-character limit during `pip install`, and
OneDrive sync and `.git` don't mix well.

```powershell
git clone https://github.com/matvazp/cytosim.git $env:USERPROFILE\dev\cytosim
cd $env:USERPROFILE\dev\cytosim
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -e ".[dev]"
python -m ipykernel install --user --name cytosim --display-name "Python (cytosim)"
pytest
```

Open `notebooks/01 simple model.ipynb` in VS Code and select the
`Python (cytosim)` kernel. `.vscode/settings.json` points the workspace at
`.venv`.

## Usage

```python
from cytosim import Params, simulate, units

r = simulate(Params(q_sheath=units.ul_per_min(800), particle_diameter=units.um(10)))
print(r.f_s_required)   # Hz
```

## Units

`Params` fields are stored in SI. Type lengths as `um(...)` and flow rates as
`ul_per_s(...)` / `ml_per_min(...)` / `ul_per_min(...)`; a bare number is metres
or m³/s and will silently give nonsense.

## Instrument defaults

Taken from the current real-life instrument: 200 × 200 µm channel; sheath
20 mL/min at ~1.5 bar and ~37 °C; sample 10 µL/s diluted 1:1 with EB
(`Params.eb()`, default) or 7 µL/s diluted 1:3.2 with PR2 (`Params.pr2()`);
oval 450 nm spot 20 µm (along flow) × 100 µm, >40 mW. The observed focused
core width of 10–30 µm is used as a validation target (`tests/`), and milk
particle populations are listed in `params.MILK_PARTICLES`.

## Model assumptions (v0.2)

- Fully developed laminar flow in the rectangular duct (`flow.RectangularDuct`,
  series solution; centreline = 2.096 × mean for a square). The sample core is
  sized by mass conservation with that field (`geometry.core_size`): flow
  through the core = Q_sample. This reproduces the measured 10–30 µm core
  width with no tunable profile parameter; plug flow (38 µm) does not.
- Core cross-section: ellipse, circular by default (`core_aspect = 1`);
  `core_aspect > 1` gives a ribbon thin across the laser beam and tall along it.
- Particles across the core see velocities from the core edge to the axis
  (~1.5 % spread for the default core); the sampling rate follows from the
  fastest particle. Plug flow is reported alongside as the lower bound since
  the instrument runs at Re ≈ 2500 (transitional).
- Gaussian laser spot along the flow direction (1/e² diameter = `spot_height`).
- Particle modelled as a uniformly emitting top-hat of diameter `particle_diameter`;
  pulse = top-hat convolved with the Gaussian beam.
- Sample rate = max(`samples_per_fwhm` / FWHM, `oversampling` × 2 × f_3dB).
- Not yet modelled: detector aperture/collection optics, laser power / photon
  budget, noise, event rate and coincidence (particle concentrations are
  stored for this), transitional-flow effects (the instrument runs at Re ≈ 2500).
