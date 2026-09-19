"""Step 1: EMD -> MRC (2D images only), MRC vertically+horizontally flipped.

Run with the data folder as cwd:  python process_emd.py
Reads *.emd in cwd (non-recursive), writes MRC/<stem>.mrc (+ .mrc.txt sidecar).
No PNG/SVG here -- those come from mrc_to_scalebar.py (reads the MRC, not the EMD).

MRC origin is bottom-left, HyperSpy/numpy origin is top-left, so the MRC is written
with flip="both". mrc_to_scalebar.py undoes it (both axes; a one-axis undo would
mirror lattice images).
"""
import glob
import math
import sys
import tempfile
import traceback
from pathlib import Path

import numpy as np

# ===== CONFIG =====
INPUT_DIR = Path.cwd()
OUT_PNGSVG_DIR = INPUT_DIR / "Scalebar"
OUT_MRC_DIR = INPUT_DIR / "MRC"

CMAP = "gray"
DPI = 100

SCALEBAR_LENGTH_FRACTION = 0.25
SCALEBAR_HEIGHT_FRACTION = 0.015
SCALEBAR_MARGIN_FRACTION = 0.04
SCALEBAR_COLOR = "white"
FONT_SIZE_PER_1024PX = 25      # pt at 1024 px width; scaled by width so labels look equal across sizes
SCALEBAR_FONT_WEIGHT = "bold"
SCALEBAR_LABEL_GAP_FRACTION = 0.018

MRC_OVERWRITE = True
MRC_FLIP = "both"
INCLUDE_NON_IMAGES = "--include-spectra" in sys.argv   # EDS cubes / spectra are skipped by default

_MRC_OK = {"int8", "int16", "uint8", "uint16", "float32", "complex64"}
_UNIT_TO_NM = {
    "nm": 1.0, "nanometer": 1.0, "nanometre": 1.0,
    "µm": 1e3, "um": 1e3, "micron": 1e3, "micrometer": 1e3, "micrometre": 1e3,
    "Å": 0.1, "a": 0.1, "angstrom": 0.1, "ångström": 0.1,
    "pm": 1e-3, "mm": 1e6,
}


# ===== scale bar helpers (imported by mrc_to_scalebar.py) =====
def _nice(x):
    e = math.floor(math.log10(x))
    b = x / 10 ** e
    n = 1.0 if b < 2 else (2.0 if b < 5 else 5.0)
    return n * 10 ** e


def _font_size_pt(w):
    return FONT_SIZE_PER_1024PX * w / 1024


def _scalebar(pixel_size, n_px, units):
    raw = pixel_size * n_px * SCALEBAR_LENGTH_FRACTION
    factor = _UNIT_TO_NM.get(str(units).strip().lower())
    if factor is None:
        val = _nice(raw)
        return val / pixel_size, f"{val:g} {units}"
    nm = _nice(raw * factor)
    bar_px = (nm / factor) / pixel_size
    label = f"{nm / 1000:g} µm" if nm >= 1000 else f"{nm:g} nm"
    return bar_px, label


# ===== MRC helpers =====
def _safe(name):
    out = "".join(c if c.isalnum() or c in "-_." else "_" for c in name.strip())
    return out or "signal"


def _voxel_size(sig):
    am = sig.axes_manager
    sig_ax = list(am.signal_axes)
    nav_ax = list(am.navigation_axes)
    if sig.data.ndim == 2 and len(sig_ax) == 2:
        sx, sy, sz = sig_ax[0].scale, sig_ax[1].scale, 1.0
    elif sig.data.ndim == 3:
        sx = nav_ax[0].scale if len(nav_ax) > 0 else 1.0
        sy = nav_ax[1].scale if len(nav_ax) > 1 else 1.0
        sz = sig_ax[0].scale if len(sig_ax) == 1 else 1.0
    else:
        sx = sy = sz = 1.0
    return tuple(float(v) if np.isfinite(v) and v > 0 else 1.0 for v in (sx, sy, sz))


def _spatial_flip(arr, mode):
    if mode == "none" or arr.ndim < 2:
        return arr
    if mode in ("both", "ud"):
        arr = arr[::-1, ...]
    if mode in ("both", "lr"):
        arr = arr[:, ::-1, ...]
    return np.ascontiguousarray(arr)


def _sidecar_text(sig, src, arr, flip):
    lines = [
        f"source          : {src}",
        f"signal title    : {sig.metadata.General.get_item('title', '')}",
        f"signal type     : {sig.metadata.Signal.get_item('signal_type', '')}",
        f"array shape     : {arr.shape}",
        f"dtype in -> out : {sig.data.dtype} -> {arr.dtype}",
        f"spatial flip    : {flip}",
        "axes:",
    ]
    for ax in sig.axes_manager._axes:
        lines.append(
            f"  {ax.name!r:>12}  size={ax.size:<6} scale={ax.scale:<12g} "
            f"offset={ax.offset:<12g} units={ax.units!r}"
        )
    return "\n".join(lines) + "\n"


_MRC_BIG_BYTES = 2 * 1024 ** 3  # above this: mmap write (Windows write() caps at ~2GB, std() temp = RAM blowup)


def _mrc_write(out_path, arr, voxel_size, is_stack):
    import mrcfile

    if arr.nbytes <= _MRC_BIG_BYTES:
        with mrcfile.new(str(out_path), overwrite=True) as m:
            m.set_data(arr)
            if is_stack:
                m.set_image_stack()
            m.voxel_size = voxel_size
        return

    mode = mrcfile.utils.mode_from_dtype(arr.dtype)
    with mrcfile.new_mmap(str(out_path), shape=arr.shape, mrc_mode=mode, overwrite=True) as m:
        m.data[:] = arr
        if is_stack:
            m.set_image_stack()
        m.voxel_size = voxel_size
        m.reset_header_stats()


def convert_signal(sig, out_path, overwrite=False, flip="both"):
    out_path = Path(out_path)
    if out_path.exists() and not overwrite:
        return "skip"

    is_spectrum = sig.data.ndim == 1
    arr = np.ascontiguousarray(sig.data)
    if arr.dtype.name not in _MRC_OK:
        arr = arr.astype(np.float32)
    if is_spectrum:
        arr = arr.reshape(1, -1)
    else:
        arr = _spatial_flip(arr, flip)

    _mrc_write(out_path, arr, _voxel_size(sig), arr.ndim == 3)

    out_path.with_suffix(out_path.suffix + ".txt").write_text(
        _sidecar_text(sig, out_path.name, arr, "none" if is_spectrum else flip),
        encoding="utf-8",
    )
    return "ok"


def selftest():
    import hyperspy.api as hs
    import mrcfile

    data = (np.random.rand(16, 32, 48) * 1000).astype("float32")
    s = hs.signals.Signal1D(data)
    s.axes_manager[0].scale, s.axes_manager[0].units = 2.5, "nm"
    s.axes_manager[1].scale, s.axes_manager[1].units = 2.5, "nm"
    s.axes_manager[2].scale, s.axes_manager[2].units = 10.0, "eV"
    s.metadata.General.title = "selftest"

    with tempfile.TemporaryDirectory() as d:
        out = Path(d) / "st_both.mrc"
        assert convert_signal(s, out, flip="both") == "ok"
        with mrcfile.open(str(out)) as m:
            assert np.allclose(m.data, data[::-1, ::-1, :])
            vs = m.voxel_size
            assert (round(float(vs.x), 3), round(float(vs.y), 3), round(float(vs.z), 3)) == (2.5, 2.5, 10.0)
    print("selftest OK")


# ===== per-file driver =====
def process_file(path):
    import hyperspy.api as hs

    loaded = hs.load(str(path))
    sigs = loaded if isinstance(loaded, (list, tuple)) else [loaded]
    multi = len(sigs) > 1
    results = []
    for i, sig in enumerate(sigs):
        if multi:
            title = _safe(sig.metadata.General.get_item("title", f"sig{i}"))
            stem = f"{path.stem}__{i}_{title}"
        else:
            stem = path.stem

        if sig.data.ndim != 2 and not INCLUDE_NON_IMAGES:
            results.append((stem, "skip (not a 2D image)"))
            continue
        try:
            status = convert_signal(sig, OUT_MRC_DIR / f"{stem}.mrc", MRC_OVERWRITE, MRC_FLIP)
            results.append((stem, f"mrc {status}"))
        except Exception as e:
            results.append((stem, f"mrc FAIL: {e}"))
            traceback.print_exc()
    return results


def main():
    OUT_MRC_DIR.mkdir(exist_ok=True)
    selftest()

    files = sorted(Path(p) for p in glob.glob(str(INPUT_DIR / "*.emd")))
    print(f"{len(files)}개 .emd 발견 ({INPUT_DIR})", flush=True)

    ok = fail = 0
    for i, fp in enumerate(files, 1):
        print(f"\n[{i}/{len(files)}] {fp.name}", flush=True)
        try:
            for stem, status in process_file(fp):
                print(f"  {status:<24} {stem}", flush=True)
                if "FAIL" in status:
                    fail += 1
                else:
                    ok += 1
        except Exception as e:
            print(f"  load error: {e}")
            traceback.print_exc()
            fail += 1

    print(f"\n완료 {ok}개, 실패 {fail}개")
    print(f"MRC: {OUT_MRC_DIR}")


if __name__ == "__main__":
    main()
