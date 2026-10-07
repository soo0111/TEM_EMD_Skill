"""Step 2: MRC -> un-flip -> (HRTEM only) rotate + crop -> percentile stretch -> scale bar -> PNG + SVG.

Run with the data folder as cwd:
    python mrc_to_scalebar.py                 # no rotation
    python mrc_to_scalebar.py --ccw 2.1       # HRTEM rotated 2.1 deg counter-clockwise (on screen)
    python mrc_to_scalebar.py --ccw -2.1      # ... clockwise
Reads MRC/*.mrc, writes Scalebar/<stem>.png/.svg and the same without bar in NoScalebar/.
MRC is un-flipped on both axes (process_emd.py flips both). HRTEM = filename contains HRTEM_MARK ("Camera Ceta").
MRC voxel size is in the EMD axis unit (nm, or µm at low mag); the unit is read from the .mrc.txt sidecar.
"""
import argparse
import glob
import os
import re
import tempfile
import traceback
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndi
from matplotlib.patches import Rectangle

from process_emd import (
    _UNIT_TO_NM, _font_size_pt, _scalebar, _spatial_flip,
    CMAP, DPI, OUT_MRC_DIR, OUT_NOBAR_DIR, OUT_PNGSVG_DIR,
    SCALEBAR_COLOR, SCALEBAR_FONT_WEIGHT,
    SCALEBAR_HEIGHT_FRACTION, SCALEBAR_LABEL_GAP_FRACTION, SCALEBAR_MARGIN_FRACTION,
)

# ===== CONFIG =====
HRTEM_MARK = "Camera Ceta"   # filename substring identifying HRTEM shots
ORDER = 3                    # bicubic, preserves lattice fringes
PCT = (1.0, 99.0)            # contrast stretch percentiles


# ===== rotation: scipy positive angle = counter-clockwise on screen =====
def largest_centred_square(rot_shape, valid, angle_deg, in_min_dim):
    rad = abs(np.radians(angle_deg))
    s = int(np.floor(in_min_dim / (np.cos(rad) + np.sin(rad)))) - 2
    cy, cx = rot_shape[0] // 2, rot_shape[1] // 2
    while s > 0:
        y0, x0 = cy - s // 2, cx - s // 2
        if valid[y0:y0 + s, x0:x0 + s].all():
            return y0, x0, s
        s -= 2
    raise RuntimeError("no valid square found")


def rotate_crop(img, angle_deg, order=ORDER):
    """Rotate by angle_deg (CCW positive), crop the largest centred square free of empty corners."""
    rot = ndi.rotate(img, angle=angle_deg, reshape=True, order=order,
                     mode="constant", cval=0.0, prefilter=True)
    mask = ndi.rotate(np.ones(img.shape, np.float32), angle=angle_deg, reshape=True,
                      order=1, mode="constant", cval=0.0, prefilter=False)
    y0, x0, s = largest_centred_square(rot.shape, mask > 0.999, angle_deg, min(img.shape))
    return rot[y0:y0 + s, x0:x0 + s].astype(np.float32), s


# ===== contrast =====
def _percentile_u8(arr, pct=PCT):
    lo, hi = np.percentile(arr, pct)
    if hi <= lo:
        hi = lo + 1.0
    return np.clip((arr - lo) / (hi - lo) * 255.0, 0, 255).astype(np.uint8)


def _is_2d_image(path):
    """Header-only peek: nz==1 => not a 3D cube; ny>1 => not a (1,N) spectrum."""
    import mrcfile

    with mrcfile.open(str(path), header_only=True, permissive=True) as m:
        return int(m.header.nz) == 1 and int(m.header.ny) > 1


def _nm_per_unit(path):
    """nm per MRC voxel-size unit, from the .mrc.txt sidecar written by process_emd.py
    (Velox gives µm at low magnification, nm at high; the MRC header has no unit)."""
    side = Path(f"{path}.txt")
    m = re.search(r"units='([^']*)'", side.read_text(encoding="utf-8")) if side.exists() else None
    unit = m.group(1).strip().lower() if m else None
    if unit in _UNIT_TO_NM:
        return _UNIT_TO_NM[unit]
    # TODO(human): no sidecar (unit is None) or unknown unit (e.g. '1/nm' diffraction) -> what to do?
    return 1.0


def load_upright(path):
    """MRC -> (float32 image with original orientation restored, pixel size nm)."""
    import mrcfile

    with mrcfile.open(str(path), mode="r", permissive=True) as m:
        arr = np.asarray(m.data, dtype=np.float32)
        px = float(m.voxel_size.x) * _nm_per_unit(path)
    return _spatial_flip(arr, "both"), px


def render_with_scalebar_arr(img_u8, pixel_size_nm, out_stem, output_dir, cmap=CMAP, bare_dir=None):
    """PNG+SVG with scale bar in output_dir; if bare_dir, also the same image without bar there."""
    h, w = img_u8.shape
    fig = plt.figure(figsize=(w / DPI, h / DPI), dpi=DPI)
    ax = plt.Axes(fig, [0.0, 0.0, 1.0, 1.0])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(img_u8, cmap=cmap, vmin=0, vmax=255, extent=(0, w, h, 0))
    if bare_dir is not None:
        for ext in ("png", "svg"):
            fig.savefig(bare_dir / f"{out_stem}.{ext}", dpi=DPI, pad_inches=0)

    bar_px, label = _scalebar(pixel_size_nm, w, "nm")
    bar_h = SCALEBAR_HEIGHT_FRACTION * h
    margin = SCALEBAR_MARGIN_FRACTION * w
    x0 = w - margin - bar_px
    y0 = h - margin - bar_h
    ax.add_patch(Rectangle((x0, y0), bar_px, bar_h,
                           facecolor=SCALEBAR_COLOR, edgecolor="none", zorder=6))
    ax.text(x0 + bar_px / 2.0, y0 - SCALEBAR_LABEL_GAP_FRACTION * h, label,
            color=SCALEBAR_COLOR, ha="center", va="bottom", zorder=6,
            fontsize=_font_size_pt(w), fontweight=SCALEBAR_FONT_WEIGHT)

    for ext in ("png", "svg"):
        fig.savefig(output_dir / f"{out_stem}.{ext}", dpi=DPI, pad_inches=0)
    plt.close(fig)


def process_one(path, ccw_deg):
    if not _is_2d_image(path):
        return "skip"
    arr, px = load_upright(path)
    if HRTEM_MARK in path.stem and ccw_deg:
        arr, _ = rotate_crop(arr, ccw_deg)
    render_with_scalebar_arr(_percentile_u8(arr), px, path.stem, OUT_PNGSVG_DIR, bare_dir=OUT_NOBAR_DIR)
    return "ok"


def selftest():
    a = np.random.rand(5, 7).astype(np.float32)
    assert np.array_equal(_spatial_flip(_spatial_flip(a, "both"), "both"), a)

    crop, s = rotate_crop(np.random.rand(64, 64).astype(np.float32), 2.1)
    assert crop.shape == (s, s) and s > 0

    # CCW positive: a bright row tilts up on the right after +angle (screen y grows downward)
    img = np.zeros((101, 101), np.float32)
    img[50, :] = 1.0
    rot = ndi.rotate(img, 10.0, reshape=False, order=1)
    assert np.argmax(rot[:, 90]) < np.argmax(rot[:, 10])

    u8 = _percentile_u8(np.random.rand(10, 10) * 1000)
    assert u8.dtype == np.uint8

    # low-mag Velox axes come in µm: pixel size must come back in nm, label in µm
    import hyperspy.api as hs
    from process_emd import convert_signal

    s = hs.signals.Signal2D(np.random.rand(32, 48).astype(np.float32))
    for a in s.axes_manager.signal_axes:
        a.scale, a.units = 0.005, "µm"
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        convert_signal(s, d / "st.mrc", overwrite=True)
        _, px = load_upright(d / "st.mrc")
        assert abs(px - 5.0) < 1e-3, px
        assert _scalebar(px, 2048, "nm")[1] == "2 µm"

        # NoScalebar copy: same image, no white bar
        for sub in ("bar", "bare"):
            (d / sub).mkdir()
        render_with_scalebar_arr(np.zeros((64, 64), np.uint8), px, "z", d / "bar", bare_dir=d / "bare")
        assert plt.imread(d / "bare" / "z.png")[..., :3].max() == 0
        assert plt.imread(d / "bar" / "z.png")[..., :3].max() == 1
        assert (d / "bare" / "z.svg").exists()
    print("selftest OK")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ccw", type=float, default=0.0, help="HRTEM rotation, degrees counter-clockwise (negative = clockwise)")
    args = ap.parse_args()

    OUT_PNGSVG_DIR.mkdir(exist_ok=True)
    OUT_NOBAR_DIR.mkdir(exist_ok=True)
    selftest()

    files = sorted(Path(p) for p in glob.glob(glob.escape(str(OUT_MRC_DIR)) + os.sep + "*.mrc"))
    print(f"{len(files)}개 .mrc 발견, HRTEM 회전: CCW {args.ccw:g} deg", flush=True)

    ok = skip = fail = 0
    for i, fp in enumerate(files, 1):
        try:
            status = process_one(fp, args.ccw)
            print(f"[{i}/{len(files)}] {status:<5} {fp.stem}", flush=True)
            if status == "ok":
                ok += 1
            else:
                skip += 1
        except Exception as e:
            print(f"[{i}/{len(files)}] FAIL  {fp.stem}: {e}")
            traceback.print_exc()
            fail += 1

    print(f"\n완료 {ok}개, 스킵 {skip}개, 실패 {fail}개")
    print(f"PNG/SVG: {OUT_PNGSVG_DIR}")
    print(f"No scale bar: {OUT_NOBAR_DIR}")


if __name__ == "__main__":
    main()
