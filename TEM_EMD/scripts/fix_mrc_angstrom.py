"""Fix MRC headers written before 2026-10-12: voxel size nm / µm -> Å (MRC standard). Pixel data untouched.

    python fix_mrc_angstrom.py <dir> [<dir> ...]     # every *.mrc under the dirs (recursive)
    python fix_mrc_angstrom.py --selftest
Only the 8 header bytes of cella.x / cella.y are rewritten, then ANGSTROM_TAG is appended to the .mrc.txt sidecar.
Idempotent: tagged sidecars are skipped, and a header that does not hold the sidecar's raw EMD scale is left alone.
"""
import re
import sys
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np

from process_emd import ANGSTROM_TAG, _UNIT_TO_NM

AXIS_RE = re.compile(r"scale=(\S+)\s+offset=\S+\s+units='([^']*)'")
CELLA_OFFSET = 40   # cella.x, cella.y, cella.z: float32 at bytes 40-51 of the 1024-byte MRC header


def fix_one(path):
    import mrcfile
    from mrcfile.utils import byte_order_from_machine_stamp

    side = Path(f"{path}.txt")
    if not side.exists():
        return "skip (no sidecar)"
    text = side.read_text(encoding="utf-8")
    if ANGSTROM_TAG in text:
        return "already done"
    axes = [(float(s), _UNIT_TO_NM[u.strip().lower()]) for s, u in AXIS_RE.findall(text)
            if u.strip().lower() in _UNIT_TO_NM]
    if not axes:
        return "skip (no length axis)"
    scale, nm_per_unit = axes[0]
    with mrcfile.open(str(path), header_only=True, permissive=True) as m:
        vx = float(m.voxel_size.x)
        cella = np.array([m.header.cella.x, m.header.cella.y], dtype=np.float64)
        bo = byte_order_from_machine_stamp(m.header.machst)
    if abs(vx / scale - 1.0) > 1e-3:
        return f"skip (header {vx:g} != sidecar scale {scale:g})"

    k = nm_per_unit * 10.0
    try:
        fh = open(path, "r+b")
    except PermissionError:
        return "skip (file in use, e.g. open in DigitalMicrograph; close it and rerun)"
    with fh:
        fh.seek(CELLA_OFFSET)
        fh.write((cella * k).astype(bo + "f4").tobytes())
    with mrcfile.open(str(path), header_only=True, permissive=True) as m:
        assert abs(float(m.voxel_size.x) / (scale * k) - 1.0) < 1e-4, f"{path}: header patch did not take"
    side.write_text(text.rstrip("\n") + "\n" + ANGSTROM_TAG + "\n", encoding="utf-8")
    return f"ok x{k:g}"


def selftest():
    import mrcfile

    with tempfile.TemporaryDirectory() as d:
        for unit, scale, want_A in (("nm", 0.0170971, 0.170971), ("µm", 0.00154921, 15.4921)):
            p = Path(d) / f"old_{unit}.mrc"
            data = np.random.rand(8, 8).astype(np.float32)
            with mrcfile.new(str(p)) as m:
                m.set_data(data)
                m.voxel_size = (scale, scale, 1.0)
            Path(f"{p}.txt").write_text(f"axes:\n           'y'  size=8   scale={scale}   offset=0   units='{unit}'\n",
                                        encoding="utf-8")
            assert fix_one(p).startswith("ok"), unit
            assert fix_one(p) == "already done", "must be idempotent"
            with mrcfile.open(str(p)) as m:
                assert abs(float(m.voxel_size.x) - want_A) < 1e-4 * want_A, (unit, float(m.voxel_size.x))
                assert np.array_equal(m.data, data), "pixel data must be untouched"
    print("selftest OK")


def main(dirs):
    selftest()
    counts = Counter()
    for d in dirs:
        for p in sorted(Path(d).rglob("*.mrc")):
            status = fix_one(p)
            counts[status.split(" (")[0].split(" x")[0]] += 1
            if not status.startswith(("ok", "already")):
                print(f"{status:<40} {p}")
    print(dict(counts))


if __name__ == "__main__":
    selftest() if sys.argv[1:] == ["--selftest"] else main(sys.argv[1:])
