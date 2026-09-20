---
name: TEM_EMD
description: Process Thermo Fisher Velox .emd files from a TEM/STEM session. Images - convert to MRC, restore orientation, save PNG + SVG with a bottom-right scale bar, rotate and crop HRTEM (Camera Ceta) images by a user-given angle after a direction preview. EDS - find the elements in each SI file, read the Velox Pre/Post-filter settings, build denoised per-element colour maps + composite + HAADF montage (from the raw SI cube if no filter was applied, or from the Velox maps as saved if the user already filtered), with the element list and colours confirmed by the user first. Use when the user has a folder of .emd files and asks for MRC/PNG/SVG export, scale bars, HRTEM rotation, or EDS colour mapping.
---

# TEM_EMD

Image pipeline: **EMD → MRC (flipped) → un-flip → [HRTEM: rotate + crop] → percentile contrast → scale bar → PNG + SVG.**
EDS pipeline: **inspect → colours (user confirms) → map** (`scripts/eds_map.py`, needs the MRC step first).
Scripts are in `scripts/` next to this file. They read the **current working directory**, so `cd` into the
data folder (the one holding the `.emd` files) and call them by absolute path.

Interpreter: the user's Python 3 (`python` or `py -3`) with
`hyperspy exspy rosettasciio h5py mrcfile scipy matplotlib numpy`
(`pip install hyperspy exspy mrcfile scipy matplotlib numpy`; check with `python -c "import hyperspy, exspy, mrcfile"` first).
Set `MPLBACKEND=Agg`. `SKILL_DIR` = the folder containing this SKILL.md
(e.g. `~/.claude/skills/TEM_EMD` or `~/.codex/skills/TEM_EMD`).

## Workflow

1. **EMD → MRC** (slow, ~3 min per 60 files — run in background):
   `python "$SKILL_DIR/scripts/process_emd.py"` → `MRC/<stem>.mrc` + `.mrc.txt`. Expect 0 FAIL.
   Multi-signal files (STEM HAADF+BF, SI) give `<stem>__<i>_<title>.mrc`; the 2-D element maps Velox saved in SI files
   (titles Zr, O, …) and the HAADF image are among them. EDS cubes/spectra are skipped
   (add `--include-spectra` only if asked). No PNG/SVG yet.
2. **HRTEM rotation — always preview first, on an image the USER picks.** HRTEM = filename contains `Camera Ceta`.
   a. Let the user choose the preview image: if they did not name one, list the `Camera Ceta` MRC file numbers
      (e.g. `ls MRC | grep "Camera Ceta" | grep -v txt`) and ask which number(s) to use (use your tool's
      choice prompt if it has one, else ask in chat; offer a few options at different magnifications).
      Never silently pick one.
   b. When the user gives a direction and angle (e.g. "CCW 2.1°"), **also render the opposite direction**:
      `python "$SKILL_DIR/scripts/preview_rotation.py" <file-number e.g. 0023> <deg>` →
      `RotationPreview/<n>_CCW<deg>.png`, `<n>_CW<deg>.png`, `<n>_compare.png` (left = CCW, right = CW).
      Run it once per chosen image.
   c. **Show the result in the chat**: display `<n>_compare.png` inline (Claude Code: `Read` the PNG; use the two
      single PNGs too if the user wants a closer look). Say which side looks right (interface / substrate edge
      horizontal, fringes level) and **wait for the user to confirm** direction and angle.
   d. If the user gave no angle, ask; if they want other angles or another image, rerun and show again.
3. **Final render** with the confirmed angle (negative = clockwise):
   `python "$SKILL_DIR/scripts/mrc_to_scalebar.py" --ccw 2.1` → `Scalebar/<stem>.png` + `.svg`.
   Omit `--ccw` when the session has no HRTEM or the user wants no rotation.
4. **Verify**: log shows 0 FAIL; ok + skip = number of `.mrc`; PNG count = SVG count = ok. Open one STEM and one
   HRTEM PNG: upright, scale bar bottom-right and readable, no black corners after rotation.

## EDS colour maps (only when the folder has `SI` .emd files and the user wants EDS maps)

Everything writes to `EDS Data/`. Run `python "$SKILL_DIR/scripts/eds_map.py" <subcommand>` (each run starts with a selftest).

5. **`inspect`** — what is in each SI file: Velox **Pre-/Post-filter** (type, size/sigma, enabled), elements Velox mapped and
   quantified, counts per pixel per element window, frames, pixel size, recommended mode, and peaks Velox did *not* map.
   Tell the user the element list and the filter state per file. Unmapped peaks are only **candidates** (an element is listed
   only if a second independent line confirms it; spectral overlaps can still fool it, and Cu/Ta/Hf/Ga are often stray
   X-rays from grid, holder, column or FIB) — ask before adding any (`--add`). Light single-line elements (B, N, F) are never
   auto-detected; add them by hand if the user wants them.
6. **`colors`** — proposes a colour per element (`EDS Data/eds_colors.json`, `confirmed: false`) and renders
   `EDS Data/color_preview.png` (swatches + sample composites + warnings). **Always show the preview inline** and ask the
   user to confirm: (a) which elements to map, (b) which go into the composite (elements added from unmapped peaks are map-only
   by default), (c) the colours. Colours follow the user's usual scheme (Hf orange, O yellow, Si green, C magenta, Zr red,
   Pt cyan); other elements get the most distinct bright colour. Relay the warnings (e.g. Zr+Si overlap looks like O; red-green
   colour blindness). Apply edits with `--elements Zr,O,Si,C`, `--set Zr=#FF0000`, `--drop Ti`, `--add Ti`,
   `--composite Zr,O,Si,C`, rerun, show again. Run **`colors --confirm` only after the user says OK**; `map` refuses to run
   without it. Do not change the colours yourself after confirmation.
7. **Filters — decide by what `inspect` found** (`map --mode auto`):
   - **Velox filter enabled and Velox maps exist → `velox` mode**: the saved maps are only coloured, never smoothed again
     (`--sigma/--post` are rejected to avoid double smoothing; `--force-extra` overrides). Use this when the user "already applied
     filters". Elements without a Velox map are skipped with a warning.
   - **No Velox filter (raw), or `--mode raw` → `raw` mode**: the map is rebuilt from the SI cube in the EMD. Pre = Gaussian
     (sigma chosen from the counts so the median element gets ~8 counts under the kernel, capped at 0.5 nm FWHM, same sigma for all files
     of a run), Post = Average 3. The filter is the agent's call; optionally show `filters [--stem 0058]` →
     `EDS Data/filter_preview.png` (Velox map vs sigma candidates with FWHM and expected noise) and say which sigma you use.
     Override with `--sigma`, `--post`, `--frames a:b` (Velox frame range), `--maps intensity`.
8. **`map`** → `EDS Data/<Velox_Pre…_Post… | Raw_PreGauss…_PostAvg…>/`: per element PNG+SVG, Additive and MaxProjection composites,
   `<stem>__Montage` (HAADF | each element | composite, paper-style, PNG+SVG) and `<stem>__report.txt/json`
   (mode, filters, counts/px, expected noise, calibration r, warnings, and the **equivalent Velox settings** so the user can re-export
   from Velox). Open the montage and one composite, read the warnings, and tell the user: mode, filter values, warnings.
9. `map` needs the MRC folder for Velox maps, HAADF and weight calibration; raw mode still runs without it (uncalibrated, no HAADF panel).

## Rules baked into the scripts (don't change without being asked)

- MRC origin is bottom-left, so `process_emd.py` writes it flipped on **both** axes (up-down and left-right).
  Restoring must also flip **both** axes; flipping only up-down mirrors lattice images. The SI cube read straight from the EMD is
  already upright.
- PNG/SVG always come from the **MRC**, never straight from EMD (direct EMD render was low-contrast/blurry).
- Contrast: 1–99 percentile stretch. Scale bar: white bar, bold label, bottom-right, length ≈ 25% width
  rounded to 1/2/5×10ⁿ nm, font size ∝ image width (25 pt at 1024 px) so labels look the same across sizes.
- Rotation: bicubic; scipy positive angle = counter-clockwise on screen; then the largest centred square without
  empty corners is cropped (so rotated HRTEM outputs are square and slightly smaller; the scale bar is
  recomputed for the cropped width). STEM images are never rotated.
- Scripts must run with cwd = data folder; outputs go to `MRC/`, `Scalebar/`, `RotationPreview/`, `EDS Data/` there.
- EDS: Velox's element maps are **composition fractions** (per pixel the quantified elements sum to 1), so raw mode builds
  fractions too (weights fitted to the Velox maps of the same file, else 1 = "relative-intensity fractions, uncalibrated").
  Velox filters are spatial blurs: Pre (on the per-pixel spectra, before quantification) and Post (on the maps); UI types
  Average / Gaussian (+ Radial Wiener for Post). They live in the EMD under `SharedProperties/SpectralFilterSettings` and
  `UnaryImageFilterSettings` (read with h5py; hyperspy does not expose them).

## Gotchas

- 2-D elemental maps inside SI files (Zr, O, …) also become grey MRC/PNG in step 3; the colour versions come from `eds_map.py`.
- Long MRC step: log is buffered — check `MRC/` file count for progress.
- Mixed pixel sizes are fine: each image uses its own MRC voxel size for the scale bar.
- Assumes Thermo Fisher Velox `.emd` and `Camera Ceta` in HRTEM filenames; change `HRTEM_MARK` in `mrc_to_scalebar.py` for other cameras.
- The rotation angle is session-specific (one session needed CCW 2.1°) — always preview, never reuse blindly.
- **EDS maps are qualitative.** Counts are very sparse (often 0.02–0.1 per pixel per element), so even filtered maps show ~30–60 %
  counting noise; C (and any element below ~0.03 counts/px) stays noisy. Raw mode reproduces the Velox maps well for Si/C/O
  (r ≈ 0.9–0.98) but not for elements whose lines overlap others (Zr r ≈ 0.6, Pt next to Zr much worse); the report warns when
  r < 0.7. Absolute fractions shift with the filter strength, so do not compare them with Velox numbers.
- Raw mode loads the whole SI cube (~1 GB RAM for a 300×500 map) and takes ~20 s per file; `--frames` drops drifting frames.
- The auto sigma is capped by a 0.5 nm blur budget; at low magnification (large pixels) it can be too small — pass `--sigma`.
