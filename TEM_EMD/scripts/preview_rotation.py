"""Rotation-direction preview for one HRTEM MRC: CCW +DEG vs CW DEG, side by side.

Run with the data folder as cwd (after process_emd.py):
    python preview_rotation.py 0023 2.1
Writes RotationPreview/<prefix>_CCW<deg>.png/.svg, _CW<deg>.png/.svg and <prefix>_compare.png.
"""
import glob
import sys
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt

from mrc_to_scalebar import _percentile_u8, load_upright, render_with_scalebar_arr, rotate_crop
from process_emd import INPUT_DIR, OUT_MRC_DIR

OUT_DIR = INPUT_DIR / "RotationPreview"


def main(prefix, deg):
    OUT_DIR.mkdir(exist_ok=True)
    (path,) = glob.glob(str(OUT_MRC_DIR / f"{prefix} *.mrc"))
    arr, px = load_upright(Path(path))

    tag = f"{deg:g}".replace(".", "_").replace("-", "")
    variants = {f"CCW{tag}": +deg, f"CW{tag}": -deg}
    for name, ang in variants.items():
        crop, _ = rotate_crop(arr, ang)
        render_with_scalebar_arr(_percentile_u8(crop), px, f"{prefix}_{name}", OUT_DIR)

    fig, axs = plt.subplots(1, 2, figsize=(14, 7.4))
    for ax, (name, ang) in zip(axs, variants.items()):
        ax.imshow(mpimg.imread(OUT_DIR / f"{prefix}_{name}.png"))
        ax.set_title(f"{'CCW' if ang > 0 else 'CW'} {abs(deg):g}°")
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(OUT_DIR / f"{prefix}_compare.png", dpi=100)
    print("saved", OUT_DIR / f"{prefix}_compare.png")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: python preview_rotation.py <file-number-prefix e.g. 0023> <degrees e.g. 2.1>")
    main(sys.argv[1], float(sys.argv[2]))
