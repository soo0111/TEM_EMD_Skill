"""EDS colour maps from Velox .emd files (elements, filters, colours, montage).

Run with the data folder as cwd (after process_emd.py, which gives MRC/ with the Velox maps + HAADF):
    python eds_map.py inspect                       # elements, Velox filters, counts per SI file -> EDS Data/eds_inspect.json
    python eds_map.py colors [--set Zr=#FF0000] [--drop Ti] [--add Ti] [--composite Zr,O,Si,C]   # draft + preview
    python eds_map.py colors --confirm              # only after the user agreed to elements + colours
    python eds_map.py filters [--stem 0058]         # raw mode: compare sigma candidates side by side
    python eds_map.py map [--mode auto|velox|raw] [--sigma S] [--post N] [--frames a:b] [--maps fraction|intensity]

Velox pipeline: Pre-filter (spatial blur of the counts) -> quantification -> Post-filter (blur of the maps).
  velox mode: the element maps Velox saved (already filtered) are only coloured, never smoothed again.
  raw mode  : rebuilt from the SI cube in the EMD: window sums -> Pre Gaussian -> fractions -> Post Average.
              Fraction weights are fitted to the Velox maps of the same file when they exist (else 1 = uncalibrated).
Maps are qualitative, not quantified compositions.
"""
import argparse
import colorsys
import glob
import json
import math
import sys
import traceback
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, to_hex, to_rgb  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from scipy.ndimage import gaussian_filter, uniform_filter  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402

from mrc_to_scalebar import _percentile_u8, load_upright, render_with_scalebar_arr  # noqa: E402
from process_emd import (  # noqa: E402
    DPI, INPUT_DIR, OUT_MRC_DIR, SCALEBAR_COLOR, SCALEBAR_FONT_WEIGHT, SCALEBAR_HEIGHT_FRACTION,
    SCALEBAR_LABEL_GAP_FRACTION, SCALEBAR_MARGIN_FRACTION, _UNIT_TO_NM, _font_size_pt, _scalebar,
)

# ===== CONFIG =====
OUT_DIR = INPUT_DIR / "EDS Data"
INSPECT_JSON = OUT_DIR / "eds_inspect.json"
COLORS_JSON = OUT_DIR / "eds_colors.json"
TARGET_COUNTS = 8        # counts under the Pre-filter kernel wanted for the median element (noise ~35 %)
MAX_FWHM_NM = 0.5        # total blur budget; interfaces of interest are ~1-3 nm wide
POST_SIZE = 3            # Post-filter Average size (Velox: odd only)
LINE_RANGE_KEV = (0.15, 12.0)
DETECT_Z, DETECT_MIN_COUNTS = 6.0, 200
# colours the user already uses in figures; anything else gets the most distinct free colour
CONVENTION = {"Hf": "#FF8000", "O": "#FFFF00", "Si": "#00CC00", "C": "#FF00FF", "Zr": "#FF0000", "Pt": "#00FFFF"}
UI_NAME = {"TopHat": "Average", "Mean": "Average", "Gaussian": "Gaussian", "Blur": "Gaussian"}  # Blur = inferred
_MACHADO_DEUTAN = np.array([[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413],
                            [-0.011820, 0.042940, 0.968881]])
_TABLE = {}


# ===== element table (exspy) =====
def table():
    """symbol -> (Z, {line: (keV, weight)})"""
    if not _TABLE:
        import exspy.material as M

        for sym, v in M.elements.as_dictionary().items():
            if not isinstance(v, dict) or "General_properties" not in v:   # skip non-element bookkeeping keys
                continue
            lines = {n: (p["energy (keV)"], p["weight"]) for n, p in v["Atomic_properties"].get("Xray_lines", {}).items()}
            _TABLE[sym] = (v["General_properties"]["Z"], lines)
    return _TABLE


def fwhm_kev(e, mnka=130.0):
    return max(math.sqrt(max(mnka ** 2 + 2.68 * (e * 1000 - 5898), 0)) / 1000, 0.06)


def main_lines(sym):
    return [n for n, (e, _) in table()[sym][1].items()
            if n in ("Ka", "La", "Ma") and LINE_RANGE_KEV[0] <= e <= LINE_RANGE_KEV[1]]


def line_window(sym, line, mnka):
    """keV window of the main line plus same-family neighbours with weight >= 0.3 (e.g. Zr La + Lb1)."""
    lines = table()[sym][1]
    e0 = lines[line][0]
    fam = [e0] + [e for n, (e, w) in lines.items() if n[0] == line[0] and w >= 0.3 and abs(e - e0) <= 0.25]
    lo, hi = min(fam), max(fam)
    return lo - 0.75 * fwhm_kev(lo, mnka), hi + 0.75 * fwhm_kev(hi, mnka)


def _idx(ax, e):
    return int(np.clip(round((e - ax.offset) / ax.scale), 0, ax.size - 1))


def _net(spec, a, b):
    """window counts, flank-background counts, significance z from a 1-D spectrum."""
    w = b - a
    win = spec[a:b].sum()
    bkg = 0.5 * (spec[max(a - w, 0):a].sum() + spec[b:b + w].sum())
    return win, bkg, (win - bkg) / math.sqrt(max(win + 0.25 * bkg * 2, 1.0))


def best_window(spec, ax, sym, mnka):
    """strongest of Ka/La/Ma for this element in this spectrum -> (line, lo, hi, a, b, net, z) or None."""
    best = None
    for ln in main_lines(sym):
        lo, hi = line_window(sym, ln, mnka)
        a, b = _idx(ax, lo), _idx(ax, hi) + 1
        win, bkg, z = _net(spec, a, b)
        if best is None or win - bkg > best[5]:
            best = (ln, lo, hi, a, b, win - bkg, z)
    return best


def _lines(sym, min_w=0.05, emax=18.0):
    return [(n, e) for n, (e, w) in table()[sym][1].items() if w >= min_w and LINE_RANGE_KEV[0] <= e <= emax]


def detect_peaks(spec, ax, mnka, exclude=()):
    """elements not in `exclude` whose main line is significant AND confirmed by a second independent line.
    Light elements with a single line (B, N, F, ...) cannot be confirmed this way: add them by hand."""
    hits = []
    for sym in table():
        b = best_window(spec, ax, sym, mnka) if main_lines(sym) and sym not in exclude else None
        if b is None:
            continue
        e = table()[sym][1][b[0]][0]
        if b[6] >= DETECT_Z and b[5] >= DETECT_MIN_COUNTS:
            hits.append({"element": sym, "line": b[0], "energy_keV": round(e, 3), "net_counts": int(b[5]), "z": round(b[6], 1)})
    # a strong known peak leaks into neighbouring windows (Si Ka -> "Ta/Rb/W", O Ka -> "Cr/V", Zr La -> "Pt/P/Au"), so keep an
    # element only if a second, independent line of it (not within one FWHM of its main line or of any known element's line) shows up too
    known_e = [e for s in exclude for _, e in _lines(s)]
    kept = []
    for h in hits:
        for n, e in _lines(h["element"]):
            if abs(e - h["energy_keV"]) <= 2 * fwhm_kev(h["energy_keV"], mnka) or any(abs(e - k) <= fwhm_kev(e, mnka) for k in known_e):
                continue
            lo, hi = line_window(h["element"], n, mnka)
            win, bkg, z = _net(spec, _idx(ax, lo), _idx(ax, hi) + 1)
            if z >= 3 and win - bkg >= 30:
                kept.append(dict(h, confirmed_by=f"{n} {e:.2f} keV"))
                break
    return sorted(kept, key=lambda h: -h["z"])[:15]


# ===== Velox settings (h5py; hyperspy does not expose them) =====
def _first_json(sp, name):
    if name not in sp or not len(sp[name]):
        return {}
    v = sp[name][next(iter(sp[name]))][0]
    return json.loads(v.decode("utf-8") if isinstance(v, bytes) else str(v))


def _ui_filter(d, kind_key, size_key):
    if not d:
        return {"enabled": False, "label": "n/a"}
    kind = UI_NAME.get(d.get(kind_key), str(d.get(kind_key)))
    f = {"enabled": bool(d.get("enabled")), "kind": kind, "size": d.get(size_key) or 0, "sigma": d.get("sigma") or 0,
         "velox_name": d.get(kind_key)}
    f["label"] = "off" if not f["enabled"] else (f"Average {int(f['size'])}" if kind == "Average"
                                                 else f"Gaussian sigma {f['sigma']:g}" if kind == "Gaussian" else kind)
    return f


def velox_info(path):
    """None if the EMD has no spectrum image, else Pre/Post filter + quantified element symbols."""
    import h5py

    with h5py.File(path, "r") as f:
        if "Data" not in f or "SpectrumImage" not in f["Data"] or not len(f["Data"]["SpectrumImage"]):
            return None
        sp = f["SharedProperties"]
        q = _first_json(sp, "EDSSpectrumQuantificationSettings")
        pre = _ui_filter(_first_json(sp, "SpectralFilterSettings"), "type", "size")
        post = _ui_filter(_first_json(sp, "UnaryImageFilterSettings"), "filter", "kernelSize")
    z2s = {z: s for s, (z, _) in table().items()}
    return {"pre": pre, "post": post, "quantified": [z2s[z] for z in q.get("elementSelection", []) if z in z2s]}


def emd_files(stems=()):
    out = []
    for p in sorted(Path(x) for x in glob.glob(str(INPUT_DIR / "*.emd"))):
        if stems and not any(s in p.stem for s in stems):
            continue
        info = velox_info(p)
        if info:
            out.append((p, info))
    return out


def velox_maps(stem):
    """{element: MRC path} of the maps Velox saved (titles that are element symbols) + HAADF path if any."""
    maps, haadf = {}, None
    for p in sorted(Path(f) for f in glob.glob(str(OUT_MRC_DIR / (glob.escape(stem) + "__*.mrc")))):
        tok = p.stem.rsplit("_", 1)[-1]
        if tok in table():
            maps[tok] = p
        elif tok == "HAADF":
            haadf = p
    return maps, haadf


def pick_mode(info, maps):
    return "velox" if (info["pre"]["enabled"] or info["post"]["enabled"]) and maps else "raw"


# ===== raw SI cube =====
def load_cube(emd, frames=None):
    import hyperspy.api as hs

    kw = {"first_frame": frames[0], "last_frame": frames[1]} if frames else {}
    s = hs.load(str(emd), select_type="spectrum_image", lazy=True, **kw)
    s = s[0] if isinstance(s, (list, tuple)) else s
    ax = s.axes_manager.signal_axes[0]
    nx = s.axes_manager["x"]
    px = float(nx.scale) * _UNIT_TO_NM.get(str(nx.units).strip().lower(), 1.0)
    mnka = float(s.metadata.get_item("Acquisition_instrument.TEM.Detector.EDS.energy_resolution_MnKa", 130.0))
    frames_n = s.metadata.get_item("Acquisition_instrument.TEM.Detector.EDS.number_of_frames", None)
    spec = np.asarray(s.data.sum(axis=(0, 1)), dtype=np.float64)
    return {"data": s.data, "ax": ax, "px_nm": px, "mnka": mnka, "spec": spec, "n_frames": frames_n}


def window_counts(cube, syms, bkg=False):
    """{sym: (float32 map, window info)}; bkg=True subtracts the mean of the two flank windows."""
    out = {}
    for sym in syms:
        b = best_window(cube["spec"], cube["ax"], sym, cube["mnka"])
        if b is None or b[5] <= 0:
            continue
        ln, lo, hi, a, e, net, z = b
        m = np.asarray(cube["data"][..., a:e], dtype=np.float32).sum(-1)
        if bkg:
            w = e - a
            m -= 0.5 * (np.asarray(cube["data"][..., max(a - w, 0):a], np.float32).sum(-1)
                        + np.asarray(cube["data"][..., e:e + w], np.float32).sum(-1))
        out[sym] = (m, {"line": ln, "keV": [round(lo, 3), round(hi, 3)], "net": int(net), "z": round(float(z), 1)})
    return out


def home_lambda(win):
    """mean counts/pixel where the element lives (pixels in the upper half of its blurred map)."""
    sm = gaussian_filter(win, 6)
    return float(win[sm >= np.median(sm)].mean())


def sigma_auto(lams, px_nm):
    lam = float(np.median(lams))
    s = math.sqrt(TARGET_COUNTS / (4 * math.pi * lam))
    s = min(s, MAX_FWHM_NM / (2.3548 * px_nm))
    return max(round(s * 2) / 2, 1.0)


def kernel_stats(sigma, post):
    """2-D effective pixels averaged and FWHM (px) of Gaussian(sigma) followed by Average(post)."""
    r = int(math.ceil(4 * sigma))
    x = np.arange(-r, r + 1)
    g = np.exp(-x ** 2 / (2 * sigma ** 2))
    k = np.convolve(g / g.sum(), np.ones(post) / post)
    xs = np.arange(len(k)) - len(k) // 2
    xf = np.linspace(xs[0] - 1, xs[-1] + 1, 20001)
    yf = np.interp(xf, xs, k, left=0, right=0)
    above = xf[yf >= k.max() / 2]
    return (1 / np.sum(k ** 2)) ** 2, float(above[-1] - above[0])


def make_filter(f):
    if not f.get("enabled"):
        return lambda a: a
    if f["kind"] == "Average":
        return lambda a, n=int(f["size"]): uniform_filter(a, n)
    if f["kind"] == "Gaussian":
        return lambda a, s=float(f["sigma"]): gaussian_filter(a, s)
    return None  # e.g. Radial Wiener: cannot be emulated


def fractions(counts, w, pre, post):
    """counts {sym: map} -> Pre blur -> weighted fractions (sum 1 where any signal) -> Post blur."""
    v = {s: w.get(s, 1.0) * np.clip(pre(c), 0, None) for s, c in counts.items()}
    tot = sum(v.values())
    return {s: post(np.divide(x, tot, out=np.zeros_like(tot), where=tot > 0)) for s, x in v.items()}


def fit_weights(counts, saved, info):
    """weights so that emulating the file's own Velox filters reproduces the Velox maps; returns (w, r per element)."""
    pre, post = make_filter(info["pre"]), make_filter(info["post"])
    syms = [s for s in counts if s in saved]
    if pre is None or post is None or len(syms) < 2:
        return {}, {}
    sub = {s: counts[s] for s in syms}
    ref = syms[-1]

    def res(lw):
        w = dict(zip(syms[:-1], np.exp(lw)), **{ref: 1.0})
        p = fractions(sub, w, pre, post)
        return np.concatenate([(p[s] - saved[s]).ravel() for s in syms])

    lw = least_squares(res, np.zeros(len(syms) - 1), diff_step=1e-2, max_nfev=25).x
    w = dict(zip(syms[:-1], np.exp(lw)), **{ref: 1.0})
    p = fractions(sub, w, pre, post)
    return w, {s: float(np.corrcoef(p[s].ravel(), saved[s].ravel())[0, 1]) for s in syms}


# ===== colours =====
def _lab_lin(c):
    c = np.asarray(c, float)
    xyz = (c @ np.array([[0.4124, 0.2126, 0.0193], [0.3576, 0.7152, 0.1192], [0.1805, 0.0722, 0.9505]])
           / np.array([0.95047, 1.0, 1.08883]))
    f = np.where(xyz > 0.008856, np.cbrt(xyz), 7.787 * xyz + 16 / 116)
    return np.array([116 * f[1] - 16, 500 * (f[0] - f[1]), 200 * (f[1] - f[2])])


def _lin(rgb):
    c = np.asarray(rgb, float)
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def delta_e(a, b, cvd=False):
    la, lb = _lin(to_rgb(a) if isinstance(a, str) else a), _lin(to_rgb(b) if isinstance(b, str) else b)
    if cvd:
        la, lb = np.clip(_MACHADO_DEUTAN @ la, 0, 1), np.clip(_MACHADO_DEUTAN @ lb, 0, 1)
    return float(np.linalg.norm(_lab_lin(la) - _lab_lin(lb)))


def assign_colors(elements, fixed=None):
    """convention colours first, then for each remaining element the visible candidate farthest from all used colours."""
    cols = dict(fixed or {})
    for e in elements:
        if e not in cols and e in CONVENTION:
            cols[e] = CONVENTION[e]
    cands = [to_hex(colorsys.hsv_to_rgb(h / 360, s, v)) for h in range(0, 360, 15)
             for s, v in ((1, 1), (0.5, 1))]   # bright only: dark hues vanish on the black map background
    cands = [c for c in cands if _lab_lin(_lin(to_rgb(c)))[0] >= 50 and delta_e(c, "#FFFFFF") >= 30]
    for e in elements:
        if e not in cols:
            used = list(cols.values()) or ["#000000"]
            cols[e] = max(cands, key=lambda c: min(delta_e(c, u) for u in used))
    return {e: cols[e] for e in elements}


def color_warnings(cols, mix_thr=30.0, cvd_thr=15.0):
    warns, els = [], list(cols)
    for i, a in enumerate(els):
        for b in els[i + 1:]:
            mix = np.clip(np.array(to_rgb(cols[a])) + np.array(to_rgb(cols[b])), 0, 1)
            for c in els:
                if c not in (a, b) and delta_e(mix, cols[c]) < mix_thr:
                    warns.append(f"{a}+{b} overlap looks like {c} (dE {delta_e(mix, cols[c]):.0f})")
            if delta_e(cols[a], cols[b], cvd=True) < cvd_thr:
                warns.append(f"{a} and {b} are hard to tell apart with red-green colour blindness")
    return warns


def load_colors():
    return json.loads(COLORS_JSON.read_text(encoding="utf-8")) if COLORS_JSON.exists() else None


def require_confirmed(cfg):
    if not cfg or not cfg.get("confirmed"):
        sys.exit("원소/색상이 아직 확정되지 않았습니다. `colors`로 미리보기를 만들고, 사용자 확인 후 `colors --confirm` 을 실행하세요.")
    return cfg


# ===== rendering =====
def _cmap(color):
    return LinearSegmentedColormap.from_list("k_" + color, ["black", color], N=256)


def _draw_scalebar(ax, w, h, px_nm):
    bar_px, label = _scalebar(px_nm, w, "nm")
    bar_h, margin = SCALEBAR_HEIGHT_FRACTION * h, SCALEBAR_MARGIN_FRACTION * w
    x0, y0 = w - margin - bar_px, h - margin - bar_h
    ax.add_patch(Rectangle((x0, y0), bar_px, bar_h, facecolor=SCALEBAR_COLOR, edgecolor="none", zorder=6))
    ax.text(x0 + bar_px / 2, y0 - SCALEBAR_LABEL_GAP_FRACTION * h, label, color=SCALEBAR_COLOR, ha="center",
            va="bottom", zorder=6, fontsize=_font_size_pt(w), fontweight=SCALEBAR_FONT_WEIGHT)


def render_rgb(rgb, px_nm, out_stem, out_dir):
    h, w = rgb.shape[:2]
    fig = plt.figure(figsize=(w / DPI, h / DPI), dpi=DPI)
    ax = plt.Axes(fig, [0, 0, 1, 1])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(rgb, extent=(0, w, h, 0))
    _draw_scalebar(ax, w, h, px_nm)
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"{out_stem}.{ext}", dpi=DPI, pad_inches=0)
    plt.close(fig)


def composite(norm, cols, elements, mode="sum"):
    st = np.stack([norm[e][..., None] * np.array(to_rgb(cols[e]), np.float32) for e in elements])
    return np.clip(st.sum(0), 0, 1) if mode == "sum" else st.max(0)


def montage(stem, px_nm, haadf_u8, norm, cols, comp_els, out_dir):
    """paper-style sheet: HAADF | each element | composite, labels in element colours, scale bar on HAADF + composite."""
    panels = ([("HAADF", haadf_u8, "white", "gray")] if haadf_u8 is not None else [])
    panels += [(e, norm[e], cols[e], _cmap(cols[e])) for e in norm]
    panels.append(("+".join(comp_els), composite(norm, cols, comp_els), "white", None))
    n = len(panels)
    ncols = 3 if n <= 6 else 4
    nrows = math.ceil(n / ncols)
    h, w = panels[0][1].shape[:2]
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * w / DPI, nrows * h / DPI), dpi=DPI, facecolor="black")
    fig.subplots_adjust(0, 0, 1, 1, wspace=0.01, hspace=0.01)
    for ax in np.atleast_1d(axes).ravel():
        ax.axis("off")
    for ax, (label, img, lc, cm) in zip(np.atleast_1d(axes).ravel(), panels):
        if img.ndim == 2:
            ax.imshow(img, cmap=cm, vmin=0, vmax=255 if img.dtype == np.uint8 else 1, extent=(0, w, h, 0))
        else:
            ax.imshow(img, extent=(0, w, h, 0))
        ax.text(0.03 * w, 0.03 * h, label, color=lc, ha="left", va="top", fontsize=_font_size_pt(w), fontweight="bold",
                path_effects=[pe.withStroke(linewidth=3, foreground="black")])   # readable on its own bright colour
        if label == "HAADF" or img.ndim == 3:
            _draw_scalebar(ax, w, h, px_nm)
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"{stem}__Montage.{ext}", dpi=DPI, facecolor="black")
    plt.close(fig)


# ===== dataset -> maps =====
def build_maps(emd, info, mode, elements, args):
    """-> dict(maps {sym: upright float map}, px_nm, rep {...report...}); mode velox|raw."""
    stem = emd.stem
    vmaps, haadf = velox_maps(stem)
    rep = {"dataset": stem, "mode": mode, "velox_pre": info["pre"]["label"], "velox_post": info["post"]["label"],
           "warnings": []}
    if mode == "velox":
        if (args.sigma is not None or args.post is not None) and not args.force_extra:
            sys.exit("이 EMD는 Velox 필터가 이미 적용돼 있습니다. --sigma/--post를 주면 이중 스무딩이 됩니다 "
                     "(그래도 하려면 --force-extra, 처음부터 다시 하려면 --mode raw).")
        if not vmaps:
            raise RuntimeError("no Velox element maps in MRC/ for this file: run process_emd.py first or use --mode raw")
        miss = [e for e in elements if e not in vmaps]
        if miss:
            rep["warnings"].append(f"Velox map missing for {miss}: skipped (use --mode raw to include them)")
        maps, px = {}, None
        for e in elements:
            if e in vmaps:
                maps[e], px = load_upright(vmaps[e])
        rep.update(applied="none extra (Velox maps as saved)")
    else:
        cube = load_cube(emd, args.frames)
        px = cube["px_nm"]
        kind = args.maps
        wc = window_counts(cube, elements, bkg=(kind == "intensity"))
        for e in elements:
            if e not in wc:
                rep["warnings"].append(f"{e}: no usable X-ray line in the spectrum, skipped")
        counts = {e: m for e, (m, _) in wc.items()}
        lam = {e: home_lambda(np.clip(counts[e], 0, None)) for e in counts}
        # same filter for every dataset of a run (comparable figures): --sigma, else the first dataset's auto value
        sigma = args.sigma if args.sigma is not None else getattr(args, "sigma_shared", None) or sigma_auto(list(lam.values()), px)
        post = args.post if args.post is not None else POST_SIZE
        neff, fw = kernel_stats(sigma, post)
        pre_f, post_f = (lambda a: gaussian_filter(a, sigma)), (lambda a: uniform_filter(a, post))
        w, r = {}, {}
        if kind == "fraction" and vmaps:
            saved = {e: load_upright(p)[0] for e, p in vmaps.items() if e in counts}
            w, r = fit_weights(counts, saved, info)
        if kind == "fraction":
            maps = fractions(counts, w, pre_f, post_f)
        else:
            maps = {e: np.clip(post_f(pre_f(c)), 0, None) for e, c in counts.items()}
        cal = ("weights fitted to the Velox maps of this file" if w else "uncalibrated: weights = 1, relative-intensity fractions")
        rep.update(applied=f"Pre Gaussian sigma {sigma:g} px ({sigma * px:.3f} nm), Post Average {post}; maps={kind}; {cal}",
                   sigma_px=sigma, sigma_auto_px=sigma_auto(list(lam.values()), px), post=post,
                   fwhm_nm=round(fw * px, 3), n_eff=round(float(neff), 1),
                   windows={e: wc[e][1] for e in wc}, counts_per_px={e: round(v, 3) for e, v in lam.items()},
                   rel_noise_pct={e: round(100 / math.sqrt(max(lam[e] * neff, 1e-9))) for e in lam},
                   velox_equivalent=f"Pre-filter Gaussian blur, Sigma {sigma:g}; Post-filter Average, Size {post}",
                   calibration_r={e: round(v, 2) for e, v in r.items()}, frames=args.frames or "all",
                   n_frames=cube["n_frames"])
        rep["warnings"] += [f"calibration r for {e} is only {v:.2f}: treat {e} as low-fidelity" for e, v in r.items() if v < 0.7]
        rep["warnings"] += [f"{e}: ~{p}% relative counting noise even after filtering (qualitative only)"
                            for e, p in rep["rel_noise_pct"].items() if p > 50]
    return {"maps": {e: m for e, m in maps.items()}, "px_nm": px, "rep": rep, "haadf": haadf}


def out_label(rep):
    if rep["mode"] == "velox":
        return "Velox_Pre" + rep["velox_pre"].replace(" sigma ", "").replace(" ", "") + "_Post" + rep["velox_post"].replace(" ", "")
    return f"Raw_PreGauss{rep['sigma_px']:g}_PostAvg{rep['post']}"


def write_report(rep, out_dir):
    (out_dir / f"{rep['dataset']}__report.json").write_text(json.dumps(rep, indent=1, ensure_ascii=False), encoding="utf-8")
    lines = [f"{k}: {v}" for k, v in rep.items() if k not in ("warnings", "windows")]
    lines += [f"window {e}: {v}" for e, v in rep.get("windows", {}).items()]
    lines += [f"WARNING: {w}" for w in rep["warnings"]]
    (out_dir / f"{rep['dataset']}__report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


# ===== commands =====
def cmd_inspect(args):
    OUT_DIR.mkdir(exist_ok=True)
    files = emd_files()
    print(f"SI .emd {len(files)}개 ({INPUT_DIR})")
    result, proposed = {}, []
    for emd, info in files:
        vmaps, haadf = velox_maps(emd.stem)
        cube = load_cube(emd)
        known = list(dict.fromkeys(list(vmaps) + info["quantified"]))
        wc = window_counts(cube, known)
        lam = {e: round(home_lambda(np.clip(m, 0, None)), 3) for e, (m, _) in wc.items()}
        det = detect_peaks(cube["spec"], cube["ax"], cube["mnka"], exclude=known)
        mode = pick_mode(info, vmaps)
        sig = sigma_auto(list(lam.values()), cube["px_nm"]) if lam else None
        result[emd.stem] = {"pre": info["pre"]["label"], "post": info["post"]["label"], "velox_maps": list(vmaps),
                            "quantified_in_velox": info["quantified"], "counts_per_px": lam, "detected_unmapped": det,
                            "recommended_mode": mode, "recommended_sigma_px": sig, "px_nm": round(cube["px_nm"], 5),
                            "n_frames": cube["n_frames"], "mrc_found": bool(vmaps or haadf)}
        proposed += [e for e in (list(vmaps) or info["quantified"] or [d["element"] for d in det if d["z"] >= 15])
                     if e not in proposed]
        print(f"\n{emd.name}\n  Velox Pre: {info['pre']['label']} | Post: {info['post']['label']} -> mode {mode}"
              f"{'' if (vmaps or haadf) else '  [MRC/ not found: run process_emd.py first]'}")
        print(f"  Velox maps: {list(vmaps)} | quantified: {info['quantified']} | frames: {cube['n_frames']} | px {cube['px_nm']:.4f} nm")
        print(f"  counts/px (element window): {lam} | recommended raw sigma: {sig}")
        for d in det:
            print(f"  peak not mapped in Velox: {d['element']} {d['line']} {d['energy_keV']} keV, net {d['net_counts']}, "
                  f"z {d['z']}, confirmed by {d['confirmed_by']}")
    proposed.sort(key=lambda e: -table()[e][0])   # heavy -> light; --elements sets any other order
    INSPECT_JSON.write_text(json.dumps({"datasets": result, "proposed_elements": proposed}, indent=1), encoding="utf-8")
    print(f"\n제안 원소: {proposed}\n저장: {INSPECT_JSON}")


def _sample(args_elements, args):
    """first SI file -> (stem, maps, px, mode) with default filters, for previews."""
    emd, info = emd_files(args.stem)[0]
    vmaps, _ = velox_maps(emd.stem)
    mode = pick_mode(info, vmaps)
    ns = argparse.Namespace(sigma=None, post=None, force_extra=False, frames=None, maps="fraction")
    res = build_maps(emd, info, mode, args_elements, ns)
    return emd.stem, res


def _norm(maps):
    return {e: _percentile_u8(m).astype(np.float32) / 255.0 for e, m in maps.items()}


def cmd_colors(args):
    OUT_DIR.mkdir(exist_ok=True)
    if not INSPECT_JSON.exists():
        cmd_inspect(args)
    insp = json.loads(INSPECT_JSON.read_text(encoding="utf-8"))
    old = load_colors() or {}
    els = args.elements.split(",") if args.elements else old.get("elements") or insp["proposed_elements"]
    els = [e for e in els if e not in (args.drop or [])] + [e for e in (args.add or []) if e not in els]
    fixed = {e: c for e, c in (old.get("colors") or {}).items() if e in els}
    for kv in args.set or []:
        e, c = kv.split("=")
        fixed[e] = to_hex(to_rgb(c))
    cols = assign_colors(els, fixed)
    quant = {e for d in insp["datasets"].values() for e in d["quantified_in_velox"]} or set(els)
    prev = old.get("composite")   # first draft: Velox-quantified elements; later drafts keep the user's choice
    comp = args.composite.split(",") if args.composite else ([e for e in prev if e in els] if prev else [e for e in els if e in quant])
    warns = color_warnings({e: cols[e] for e in comp})
    stem, res = _sample(els, args)
    norm = _norm(res["maps"])
    comp = [e for e in comp if e in norm]
    cfg = {"elements": els, "colors": cols, "composite": comp, "confirmed": bool(args.confirm), "warnings": warns}
    COLORS_JSON.write_text(json.dumps(cfg, indent=1), encoding="utf-8")

    fig = plt.figure(figsize=(11, 4.2), facecolor="black")
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 3.2], hspace=0.05, wspace=0.03, left=0.02, right=0.98, top=0.97, bottom=0.03)
    axc = fig.add_subplot(gs[0, :])
    axc.set_axis_off()
    axc.set_xlim(0, len(els)); axc.set_ylim(0, 1)
    for i, e in enumerate(els):
        axc.add_patch(Rectangle((i + 0.05, 0.1), 0.9, 0.8, facecolor=cols[e], edgecolor="white", lw=0.5))
        axc.text(i + 0.5, 0.5, f"{e}{'' if e in comp else ' (map only)'}", ha="center", va="center", fontsize=11,
                 color="black" if delta_e(cols[e], "#000000") > 60 else "white", fontweight="bold")
    for j, (mode, ttl) in enumerate((("sum", "Additive"), ("max", "MaxProjection"))):
        ax = fig.add_subplot(gs[1, j]); ax.axis("off")
        ax.imshow(composite(norm, cols, comp, mode)); ax.set_title(f"{stem[:4]} {ttl}", color="white", fontsize=9)
    ax = fig.add_subplot(gs[1, 2]); ax.axis("off")
    ax.text(0, 1, "\n".join(["mode: " + res["rep"]["mode"], "composite: " + "+".join(comp), ""]
                            + (["WARNINGS"] + warns if warns else ["no colour warnings"])), color="white", va="top",
            fontsize=8, wrap=True)
    fig.savefig(OUT_DIR / "color_preview.png", dpi=110, facecolor="black")
    plt.close(fig)
    print(json.dumps(cfg, indent=1))
    print(f"미리보기: {OUT_DIR / 'color_preview.png'}")
    print("확정됨" if args.confirm else "아직 확정 아님 -> 사용자 확인 후 `colors --confirm`")


def cmd_filters(args):
    cfg = load_colors() or {}
    emd, info = emd_files(args.stem)[0]
    vmaps, _ = velox_maps(emd.stem)
    els = cfg.get("elements") or list(vmaps) or info["quantified"]
    cols = assign_colors(els, cfg.get("colors"))
    comp = [e for e in cfg.get("composite", els) if e in els]
    cube = load_cube(emd, args.frames)
    wc = window_counts(cube, els)
    counts = {e: m for e, (m, _) in wc.items()}
    lam = {e: home_lambda(np.clip(counts[e], 0, None)) for e in counts}
    s0 = sigma_auto(list(lam.values()), cube["px_nm"])
    cands = sorted({max(s0 - 1, 1.0), s0, s0 + 1})
    saved = {e: load_upright(p)[0] for e, p in vmaps.items() if e in counts}
    w, _ = fit_weights(counts, saved, info) if vmaps else ({}, {})
    panels = []
    if saved:
        panels.append((f"Velox saved (Pre {info['pre']['label']} / Post {info['post']['label']})", composite(_norm(saved), cols, [e for e in comp if e in saved])))
    for s in cands:
        neff, fw = kernel_stats(s, POST_SIZE)
        f = fractions(counts, w, lambda a, s=s: gaussian_filter(a, s), lambda a: uniform_filter(a, POST_SIZE))
        noise = " ".join(f"{e}{round(100 / math.sqrt(lam[e] * neff))}%" for e in counts)
        panels.append((f"Pre Gauss {s:g} / Post Avg {POST_SIZE}{'  <- auto' if s == s0 else ''}\nFWHM {fw * cube['px_nm']:.2f} nm  noise {noise}",
                       composite(_norm(f), cols, [e for e in comp if e in f])))
    fig, axes = plt.subplots(1, len(panels), figsize=(3.6 * len(panels), 2.9), facecolor="black")
    for ax, (t, img) in zip(np.atleast_1d(axes), panels):
        ax.imshow(img); ax.axis("off"); ax.set_title(t, color="white", fontsize=7)
    OUT_DIR.mkdir(exist_ok=True)
    fig.savefig(OUT_DIR / "filter_preview.png", dpi=110, facecolor="black", bbox_inches="tight")
    plt.close(fig)
    print(f"{emd.stem}: 자동 sigma {s0:g} px, 후보 {cands} -> {OUT_DIR / 'filter_preview.png'}")


def cmd_map(args):
    cfg = require_confirmed(load_colors())
    els, cols, comp = cfg["elements"], cfg["colors"], cfg["composite"]
    ok = fail = 0
    for emd, info in emd_files(args.stem):
        try:
            vmaps, haadf = velox_maps(emd.stem)
            mode = pick_mode(info, vmaps) if args.mode == "auto" else args.mode
            print(f"\n{emd.name}: mode {mode} (Velox Pre {info['pre']['label']}, Post {info['post']['label']})", flush=True)
            res = build_maps(emd, info, mode, els, args)
            rep, px, maps = res["rep"], res["px_nm"], res["maps"]
            if mode == "raw":
                args.sigma_shared = rep["sigma_px"]
            out = OUT_DIR / out_label(rep)
            out.mkdir(parents=True, exist_ok=True)
            norm = _norm(maps)
            for e, m in maps.items():
                render_with_scalebar_arr(_percentile_u8(m), px, f"{emd.stem}__{e}", out, cmap=_cmap(cols[e]))
            ce = [e for e in comp if e in norm]
            for m, lab in (("sum", "Additive"), ("max", "MaxProjection")):
                render_rgb(composite(norm, cols, ce, m), px, f"{emd.stem}__Composite_{''.join(ce)}_{lab}", out)
            h = None
            if res["haadf"] is not None:
                ha, _ = load_upright(res["haadf"])
                h = _percentile_u8(ha) if ha.shape == next(iter(maps.values())).shape else None
            montage(emd.stem, px, h, norm, cols, ce, out)
            write_report(rep, out)
            for w in rep["warnings"]:
                print("  WARNING:", w)
            print(f"  ok -> {out}  | {rep.get('applied')}")
            ok += 1
        except SystemExit:
            raise
        except Exception as e:
            print(f"FAIL  {emd.name}: {e}")
            traceback.print_exc()
            fail += 1
    print(f"\n완료 {ok}개 파일, 실패 {fail}개")


def selftest():
    rng = np.random.default_rng(0)
    counts = {e: rng.poisson(0.1, (40, 60)).astype(np.float32) for e in ("Zr", "O", "Si", "C")}
    assert abs(gaussian_filter(counts["Si"], 3).sum() - counts["Si"].sum()) < 1e-2 * counts["Si"].sum()  # blur keeps counts
    s = sum(fractions(counts, {}, lambda a: gaussian_filter(a, 2), lambda a: uniform_filter(a, 3)).values())
    assert np.allclose(s, 1.0, atol=1e-4)                                         # fractions sum to 1
    assert sigma_auto([0.073, 0.059, 0.118, 0.025], 0.04836) == 3.0                # rule reproduces the tuned value
    z2s = {z: s for s, (z, _) in table().items()}
    assert [z2s[z] for z in (6, 8, 14, 40)] == ["C", "O", "Si", "Zr"]
    lo, hi = line_window("Zr", "La", 130.0)
    assert lo < 2.042 < 2.124 < hi                                                 # Zr La window includes Lb1
    cols = assign_colors(["Zr", "O", "Si", "C", "Ti"])
    assert cols["Zr"] == "#FF0000" and cols["Ti"] not in CONVENTION.values()
    assert color_warnings({"Zr": "#FF0000", "Si": "#00CC00", "X": "#FFCC00"})       # red+green ~ yellow flagged
    for bad in (None, {"confirmed": False}):
        try:
            require_confirmed(bad)
            raise AssertionError("gate did not stop")
        except SystemExit:
            pass
    print("selftest OK")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("inspect")
    c = sub.add_parser("colors")
    c.add_argument("--elements"); c.add_argument("--set", action="append"); c.add_argument("--drop", action="append")
    c.add_argument("--add", action="append"); c.add_argument("--composite"); c.add_argument("--confirm", action="store_true")
    c.add_argument("--stem", action="append", default=[])
    f = sub.add_parser("filters")
    f.add_argument("--stem", action="append", default=[]); f.add_argument("--frames", type=lambda s: tuple(int(x) for x in s.split(":")))
    m = sub.add_parser("map")
    m.add_argument("--mode", choices=["auto", "velox", "raw"], default="auto")
    m.add_argument("--sigma", type=float); m.add_argument("--post", type=int)
    m.add_argument("--frames", type=lambda s: tuple(int(x) for x in s.split(":")), help="first:last EDS frame (Velox frame range)")
    m.add_argument("--maps", choices=["fraction", "intensity"], default="fraction")
    m.add_argument("--stem", action="append", default=[]); m.add_argument("--force-extra", action="store_true")
    args = ap.parse_args()
    selftest()
    {"inspect": cmd_inspect, "colors": cmd_colors, "filters": cmd_filters, "map": cmd_map}[args.cmd](args)


if __name__ == "__main__":
    main()
