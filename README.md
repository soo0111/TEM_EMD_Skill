# TEM_EMD skill

Agent skill (Claude Code / Codex) for TEM sessions saved as Thermo Fisher Velox `.emd`:

** HRTEM Images: EMD → MRC → restore orientation → (HRTEM) rotate + crop → PNG + SVG with a bottom-right scale bar.**
** STEM Images: EMD → MRC → restore orientation → PNG + SVG with a bottom-right scale bar.**

- MRC is stored flipped on both axes; the final PNG/SVG un-flips both, and are always rendered from the MRC (not straight from EMD).
- HRTEM rotation: you give a direction and angle, the skill previews **both directions** on an image you pick, you confirm, then all HRTEM images are rotated and cropped to the largest square without empty corners.
- Scale bar: white, bottom-right, label font size proportional to image width.

## Install

Copy the `TEM_EMD/` folder to:

- Claude Code: `~/.claude/skills/TEM_EMD`
- Codex: `~/.codex/skills/TEM_EMD`

Requires Python 3 with `pip install hyperspy mrcfile scipy matplotlib numpy`.

## Use

Ask the agent, e.g. 
*"Use TEM_EMD on this folder; rotate HRTEM 2° counter-clockwise."*
*"/TEM_EMD Please use this skill on this folder; rotate HRTEM 2° counter-clockwise for 23 image."*
Or run the scripts by hand from the data folder:

```bash
python TEM_EMD/scripts/process_emd.py                  # EMD -> MRC/
python TEM_EMD/scripts/preview_rotation.py 0023 2.1    # CCW vs CW preview
python TEM_EMD/scripts/mrc_to_scalebar.py --ccw 2.1    # PNG/SVG in Scalebar/ (negative = clockwise)
```

## License

MIT — see [LICENSE](LICENSE).

