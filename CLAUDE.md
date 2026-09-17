# cytosim — project context for Claude

Read this first. It records decisions and state that are not derivable from the code.

## Purpose
Simulation model to estimate the ADC sampling rate needed to digitise detector
pulses of a flow cytometer channel analysing milk (FOSS, "Pontos 2" ADC study).
Hard requirements: Python backend for the math; a frontend shareable with
colleagues; a graphical, animated view of the channel (sample core, sheath,
laser spot, detector); interactive parameters (flow rates, spot size, detector
area, particle size, ...). Sharing via SharePoint is desirable, not required.

## Decisions taken
- **Route B**: Python runs in the browser via Pyodide (stlite = Streamlit in
  the browser). Static build → host as a link (GitHub Pages / Azure Static Web
  Apps); SharePoint only links to it. Keep `src/cytosim` NumPy/SciPy-only,
  no file I/O, no plotting inside the package.
- **Option 1 dev style**: notebook (`notebooks/`) is the scratchpad for physics;
  stable formulas are promoted into `src/cytosim/` and tested. `%autoreload` on.
- Units: `Params` stores SI; inputs typed via `um()`, `ul_per_s()`, `ml_per_min()`.
  A `UserParams` in human units + `.to_si()` was discussed as a possible later
  refactor (would map 1:1 to app sliders).
- Repo: https://github.com/matvazp/cytosim (public). Working clone is
  `C:\Users\mtsv\dev\cytosim` with in-project `.venv`. The original OneDrive
  folder (`...\01 Pontos 2 ADC Sample Rate model`) is a stale copy; its `doc/`
  subfolder holds project documentation not in git.
- Git is at `C:\Program Files\Git\cmd` (may need adding to PATH in a shell).

## Model state (v0.2)
- Instrument defaults: 200×200 µm channel; sheath 20 mL/min ~37 °C; sample
  10 µL/s 1:1 EB (`Params.eb()`) or 7 µL/s 1:3.2 PR2 (`Params.pr2()`);
  450 nm oval spot 20 (along flow) × 100 µm, >40 mW.
- Velocity field: laminar rectangular-duct series solution
  (`flow.RectangularDuct`, centreline/mean = 2.096 for the square). Core sized
  by mass conservation with that field (`geometry.core_size`, elliptical core,
  `Params.core_aspect`, circular by default — injector geometry still unknown).
  `flow.PlugFlow` kept only as the transitional lower bound
  (`Results.core_diameter_plug`, `Results.f_s_required_plug`).
- Findings: Re ≈ 2500 (transitional); predicted core 26.7 µm (EB) / 22.4 µm
  (PR2), inside the measured 10–30 µm with no tunable flag (plug: 38 / 32 µm).
  The core sits on the flat top of the profile: velocity 17.7–18.0 m/s across
  it, pulse-FWHM spread only ~1.5 %, so the fastest (centreline) particle sets
  f_s. 1 µm particle → 655 ns FWHM → **~15 MHz** with the placeholder
  criterion (10 samples/FWHM, 5× Nyquist); plug-flow bound ~7 MHz. Do not
  quote until the criterion is agreed with signal processing.
- Assumed (unconfirmed): detection is fluorescence from EB-stained bacteria /
  somatic cells; fat globules are scatter/background only.
- Not modelled yet: detector aperture/collection optics, photon budget, noise,
  event rate/coincidence (`MILK_PARTICLES` concentrations are stored for it),
  transitional-flow flattening of the profile (only bracketed by plug flow).

## Next steps (in order)
1. **Event rate / coincidence** from `MILK_PARTICLES`: ~5e6 fat globules/s
   (~5 in the beam at any instant → continuous background), bacteria up to
   ~5e4/s, somatic cells ~1e4/s. Likely reframes the sampling criterion as
   "resolve two bacteria ~1 µs apart" rather than "10 samples per FWHM".
2. Agree the sampling criterion with signal processing.
3. Detector geometry + photon budget.
4. stlite app in `app/` on top of `simulate()`, with an SVG/JS animated
   channel view; deploy via CI to GitHub Pages.
Still open for the user: core shape (circular vs ribbon → `core_aspect`) and
confirmation that detection is fluorescence (EB = DNA stain).

## CI (GitHub Actions, `.github/workflows/ci.yml`)
Runs pytest on 3.11/3.12 on every push/PR. The user asked to be **proactively
reminded of CI additions when they become relevant**:
- when a second contributor appears → branch protection requiring CI on PRs;
- when code style starts to drift / contributors join → `ruff` lint step;
- when the stlite app exists → build + deploy to GitHub Pages on push to main;
- when the notebook matters as documentation → execute it in CI so it can't rot.

## Conventions
- Tests in `tests/`, run `pytest -q`. All must pass before pushing.
- Commit small, descriptive; push to `main` directly for now.
- Explain physics in plain terms when asked; the user is onboarding to the domain.
