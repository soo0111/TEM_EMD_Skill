# TEM_EMD skill

Agent skill (Claude Code / Codex) for TEM sessions saved as Thermo Fisher Velox `.emd`:

_**1. HRTEM Images: EMD → MRC → restore orientation → (HRTEM) rotate + crop → PNG + SVG with a bottom-right scale bar.**_

_**2. STEM Images: EMD → MRC → restore orientation → PNG + SVG with a bottom-right scale bar.**_

_**3. EDS (SI files): elements + Velox filter settings → confirmed colours → per-element colour maps, composite and HAADF montage.**_

- MRC is stored flipped on both axes; the final PNG/SVG un-flips both, and are always rendered from the MRC (not straight from EMD).
- HRTEM rotation: you give a direction and angle, the skill previews **both directions** on an image you pick, you confirm, then all HRTEM images are rotated and cropped to the largest square without empty corners.
- EDS: reads the Velox Pre/Post-filter settings from the EMD. Files you already filtered in Velox are coloured as saved (no double smoothing); raw files are rebuilt from the SI cube with auto-chosen filters. Elements and colours are confirmed with you before mapping. Maps are qualitative, not a quantification.
- Scale bar: white, bottom-right, label font size proportional to image width.

## Install

Copy the `TEM_EMD/` folder to:

- Claude Code: `~/.claude/skills/TEM_EMD`
- Codex: `~/.codex/skills/TEM_EMD`

Requires Python 3 with `pip install hyperspy exspy mrcfile scipy matplotlib numpy`.

## Use

Ask the agent, e.g. 

*"Use TEM_EMD on this folder; rotate HRTEM 2° counter-clockwise."*

*"/TEM_EMD Please use this skill on this folder; rotate HRTEM 2° counter-clockwise for 23 image."*

*"Use TEM_EMD to make EDS colour maps for the SI files in this folder."*

Or run the scripts by hand from the data folder:

```bash
python TEM_EMD/scripts/process_emd.py                  # EMD -> MRC/
python TEM_EMD/scripts/preview_rotation.py 0023 2.1    # CCW vs CW preview
python TEM_EMD/scripts/mrc_to_scalebar.py --ccw 2.1    # PNG/SVG in Scalebar/ (negative = clockwise)
```

EDS colour maps (after `process_emd.py`):

```bash
python TEM_EMD/scripts/eds_map.py inspect              # elements, Velox filters, counts per SI file
python TEM_EMD/scripts/eds_map.py colors               # draft colours + EDS Data/color_preview.png
python TEM_EMD/scripts/eds_map.py colors --confirm     # after you agree with elements and colours
python TEM_EMD/scripts/eds_map.py filters              # optional: compare sigma candidates (raw mode)
python TEM_EMD/scripts/eds_map.py map                  # maps, composites, montage, report in EDS Data/
```

EDS maps are qualitative (counts are sparse; not a quantification). See `TEM_EMD/SKILL.md` for the details and limits.

## License

MIT — see [LICENSE](LICENSE).

