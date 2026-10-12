# Changelog

## 2026-10-12

### Fixed
- **MRC pixel size was 10× (nm) or 10⁴× (µm) too small in other MRC readers.** The MRC standard stores the
  voxel size in Å, but `process_emd.py` wrote the raw EMD axis number (nm, or µm at low magnification).
  DigitalMicrograph, Fiji and IMOD therefore showed a 5 nm feature as 5 Å. This skill's own PNG/SVG scale
  bars were correct, because `load_upright()` read the unit from the `.mrc.txt` sidecar.
  `process_emd.py` now writes Å and marks the sidecar with `ANGSTROM_TAG`. `load_upright()` reads tagged
  files as Å and older, untagged files by their sidecar unit, so old folders still render correctly.

### Added
- `fix_mrc_angstrom.py <dir> [...]`: converts existing MRC headers to Å in place. Only the 8 bytes of `cella.x/y`
  are rewritten; pixel data is untouched. It can be rerun safely: tagged files are skipped, a header that does not
  match its sidecar is left alone, and files locked by another program (e.g. open in DigitalMicrograph) are
  listed and skipped.
- Self-checks: `process_emd.py` expects 2.5 nm → 25 Å in the header; `fix_mrc_angstrom.py` checks the nm and µm
  cases, idempotency, and that pixel data is untouched.

## 2026-10-07

### Fixed
- **Scale bar unit at low magnification.** Velox stores low-magnification pixel sizes in µm, but the MRC
  header has no unit field, and `load_upright()` always read the number as nm. Low-magnification TEM and STEM
  images (and EDS maps in `velox` mode) therefore got a bar labelled e.g. "2 nm" instead of "2 µm".
  `load_upright()` now reads the unit from the `.mrc.txt` sidecar and returns the pixel size in nm.
  Existing `MRC/` folders do not need re-conversion: rerun `mrc_to_scalebar.py` (and `eds_map.py map`).
  EDS `raw` mode was not affected (it already converted units).
- Data folders whose path contains `[` `]` were treated as glob patterns, so no `.emd` / `.mrc` files were found.

### Added
- No-scale-bar copies for editing, same file names, PNG + SVG:
  - `mrc_to_scalebar.py` → `NoScalebar/<stem>.png/.svg` next to `Scalebar/`.
  - `eds_map.py map` → `EDS Data/<run>/NoScalebar/` for the element maps and composites (not the montage).
- Self-check in `mrc_to_scalebar.py`: a µm image must come back as nm pixel size with a µm label, and the
  NoScalebar copy must contain no bar.

## 2026-09-20
- EDS colour mapping (`eds_map.py`): inspect, colours with user confirmation, Velox/raw filter modes, montage, report.

## 2026-09-19
- First version: EMD → MRC → PNG/SVG with scale bar, HRTEM rotation preview, rotate + crop.
