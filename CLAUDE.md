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

## Model state (v0.1)
- Instrument defaults: 200×200 µm channel; sheath 20 mL/min ~37 °C; sample
  10 µL/s 1:1 EB (`Params.eb()`) or 7 µL/s 1:3.2 PR2 (`Params.pr2()`);
  450 nm oval spot 20 (along flow) × 100 µm, >40 mW.
- Findings: Re ≈ 2500 (transitional); parabolic profile (core at 2× mean,
  ~17 m/s) reproduces the measured 10–30 µm core width, plug flow does not;
  1 µm particle → ~0.7 µs pulse FWHM → **~14 MHz** with the placeholder
  criterion (10 samples/FWHM, 5× Nyquist oversampling). Do not quote until
  the criterion is agreed with signal processing.
- Not modelled yet: detector aperture/collection optics, photon budget, noise,
  event rate/coincidence (`MILK_PARTICLES` concentrations are stored for it).

## Next steps (in order)
1. **Real velocity field in the square duct** (replaces the plug/parabolic
   switch, factor 2 on the answer): laminar rectangular-duct series solution
   (centreline/mean ≈ 2.10 for a square); core size/shape from mass
   conservation with that profile; velocity *distribution* across the core →
   pulse-width distribution; `f_s_required` from the fastest particle; report
   plug flow as the lower bound because Re ≈ 2500 is transitional. Test: predicted
   core width lands in the measured 10–30 µm without a tunable profile flag.
   Open questions for the user: is the core circular or a ribbon (thin across
   the laser, tall along it)? Is detection fluorescence (EB = DNA stain, so
   bacteria/somatic cells) with fat globules only as scatter/background?
2. **Event rate / coincidence** from `MILK_PARTICLES`: ~5e6 fat globules/s
   (~5 in the beam at any instant → continuous background), bacteria up to
   ~5e4/s, somatic cells ~1e4/s. Likely reframes the sampling criterion as
   "resolve two bacteria ~1 µs apart" rather than "10 samples per FWHM".
3. Agree the sampling criterion with signal processing.
4. Detector geometry + photon budget.
5. stlite app in `app/` on top of `simulate()`, with an SVG/JS animated
   channel view; deploy via CI to GitHub Pages.

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
