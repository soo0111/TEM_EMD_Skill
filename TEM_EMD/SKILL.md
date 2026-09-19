---
name: TEM_EMD
description: Convert TEM/STEM/HRTEM .emd files (Thermo Fisher Velox) to MRC, restore orientation, and save PNG + SVG with a bottom-right scale bar; rotates and crops HRTEM (Camera Ceta) images by a user-given angle after a direction preview. Use when the user has a folder of .emd files from a TEM session and asks for MRC conversion, PNG/SVG export, scale bars, flipping, or HRTEM rotation. Images only — EDS colour maps are out of scope.
---

# TEM_EMD

Pipeline: **EMD → MRC (flipped) → un-flip → [HRTEM: rotate + crop] → percentile contrast → scale bar → PNG + SVG.**
Scripts are in `scripts/` next to this file. They read the **current working directory**, so `cd` into the
data folder (the one holding the `.emd` files) and call them by absolute path.

Interpreter: the user's Python 3 (`python` or `py -3`) with `hyperspy rosettasciio mrcfile scipy matplotlib numpy`
(`pip install hyperspy mrcfile scipy matplotlib numpy`; check with `python -c "import hyperspy, mrcfile"` first).
Set `MPLBACKEND=Agg`. `SKILL_DIR` = the folder containing this SKILL.md
(e.g. `~/.claude/skills/TEM_EMD` or `~/.codex/skills/TEM_EMD`).

## Workflow

1. **EMD → MRC** (slow, ~3 min per 60 files — run in background):
   `python "$SKILL_DIR/scripts/process_emd.py"` → `MRC/<stem>.mrc` + `.mrc.txt`. Expect 0 FAIL.
   Multi-signal files (STEM HAADF+BF, SI) give `<stem>__<i>_<title>.mrc`. EDS cubes/spectra are skipped
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

## Rules baked into the scripts (don't change without being asked)

- MRC origin is bottom-left, so `process_emd.py` writes it flipped on **both** axes (up-down and left-right).
  Restoring must also flip **both** axes; flipping only up-down mirrors lattice images.
- PNG/SVG always come from the **MRC**, never straight from EMD (direct EMD render was low-contrast/blurry).
- Contrast: 1–99 percentile stretch. Scale bar: white bar, bold label, bottom-right, length ≈ 25% width
  rounded to 1/2/5×10ⁿ nm, font size ∝ image width (25 pt at 1024 px) so labels look the same across sizes.
- Rotation: bicubic; scipy positive angle = counter-clockwise on screen; then the largest centred square without
  empty corners is cropped (so rotated HRTEM outputs are square and slightly smaller; the scale bar is
  recomputed for the cropped width). STEM images are never rotated.
- Scripts must run with cwd = data folder; outputs go to `MRC/`, `Scalebar/`, `RotationPreview/` there.

## Gotchas

- 2-D elemental maps inside SI files (Zr, O, …) also become grey MRC/PNG; that is expected, colour maps are a separate task.
- Long MRC step: log is buffered — check `MRC/` file count for progress.
- Mixed pixel sizes are fine: each image uses its own MRC voxel size for the scale bar.
- Assumes Thermo Fisher Velox `.emd` and `Camera Ceta` in HRTEM filenames; change `HRTEM_MARK` in `mrc_to_scalebar.py` for other cameras.
- The rotation angle is session-specific (one session needed CCW 2.1°) — always preview, never reuse blindly.
