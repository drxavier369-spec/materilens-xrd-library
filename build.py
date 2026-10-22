"""Build a powder XRD reference library from Crystallography Open Database CIFs, in shards that can run in parallel.

    python build.py shard  cod  shard_1_0.npz  --workers 4
    python build.py merge  xrd_library_cod.npz  shard_*.npz
    python build.py plan   listing.txt  > slices.json

plan: reads `rsync --list-only -r` output of the COD cif tree and splits it into slices of similar size, each a
comma-separated list of folders such as 1/50/12 (100 COD numbers each), for the workflow's job matrix.

shard: computes the powder pattern (Cu Kα1, 5–140° 2θ, stored as d-spacings) of every inorganic CIF under the folder.
Left out:
  - COD numbers 9000000 and above, which are copied from the American Mineralogist Crystal Structure Database;
  - organic and metal-organic entries (formula has both C and H);
  - structures with more than MAX_SITES atoms in the cell, or that take longer than TIMEOUT seconds.
merge: joins shards, keeping one entry per COD number.

COD data are open and free to use (https://www.crystallography.net/cod/); each phase keeps its COD number.
Output: one .npz with flat arrays d, I, hkl (lines x 4, int16), hkl_len, offsets (per phase), elmask (2 x uint64 per
phase, bit = index in SYMBOLS) and meta (one JSON object per phase, newline-separated, starting with "id").
"""
from __future__ import annotations

import argparse
import json
import re
import signal
import warnings
from multiprocessing import Pool
from pathlib import Path

import numpy as np

AMCSD_FROM = 9000000
MAX_SITES = 400
TIMEOUT = 60
MIN_I = 0.5  # weaker lines are dropped to keep the library small

# Strongest-line intensity per unit mass of corundum, I_max / (V^2 * density) from pymatgen's unscaled pattern
# (Cu Ka1). Dividing a phase's value by this gives its RIR (I/I_corundum).
CORUNDUM_I_PER_MASS = 3.2038

SYMBOLS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y "
           "Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re "
           "Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg "
           "Cn Nh Fl Mc Lv Ts Og D").split()
Z = {s: i for i, s in enumerate(SYMBOLS)}


def reference_from_structure(structure, ref_id: str, source: str = "", licence: str = "") -> dict:
    from pymatgen.analysis.diffraction.xrd import XRDCalculator
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    calc = XRDCalculator(wavelength="CuKa1")
    pattern = calc.get_pattern(structure, two_theta_range=(5, 140))
    absolute = calc.get_pattern(structure, scaled=False, two_theta_range=(5, 140))
    rir = float(absolute.y.max()) / (structure.volume ** 2 * structure.density) / CORUNDUM_I_PER_MASS
    sga = SpacegroupAnalyzer(structure, symprec=0.1)
    lat = structure.lattice
    return {
        "id": ref_id,
        "name": structure.composition.reduced_formula,
        "formula": structure.composition.reduced_formula,
        "elements": sorted({el.symbol for el in structure.composition.elements}),
        "spacegroup": sga.get_space_group_symbol(),
        "crystal_system": sga.get_crystal_system(),
        "lattice": {k: round(v, 4) for k, v in zip("a b c alpha beta gamma".split(), lat.parameters)},
        "density": round(structure.density, 3),
        "rir": round(rir, 2),
        "source": source,
        "licence": licence,
        "d": [round(float(d), 5) for d in pattern.d_hkls],
        "I": [round(float(i), 2) for i in pattern.y],
        "hkl": [list(h[0]["hkl"]) for h in pattern.hkls],
    }


def save(refs: list[dict], path) -> None:
    meta, offsets, d, I, hkl, hl, em = [], [0], [], [], [], [], []
    for r in refs:
        meta.append(json.dumps({k: v for k, v in r.items() if k not in ("d", "I", "hkl")}, separators=(",", ":")).encode())
        offsets.append(offsets[-1] + len(r["d"]))
        d.extend(r["d"]); I.extend(r["I"])
        h = r["hkl"]
        hl.append(len(h[0]) if h else 3)
        hkl.extend(list(x) + [0] * (4 - len(x)) for x in h)
        m = np.zeros(2, np.uint64)
        for e in r["elements"]:
            if e in Z:
                m[Z[e] // 64] |= np.uint64(1) << np.uint64(Z[e] % 64)
        em.append(m)
    np.savez_compressed(path, meta=np.frombuffer(b"\n".join(meta), np.uint8), offsets=np.asarray(offsets, np.int64),
                        d=np.asarray(d, np.float32), I=np.asarray(I, np.float32),
                        hkl=np.asarray(hkl, np.int16).reshape(-1, 4), hkl_len=np.asarray(hl, np.int8),
                        elmask=np.asarray(em, np.uint64).reshape(-1, 2))


def load(path) -> list[dict]:
    z = np.load(path, allow_pickle=False)
    raw = bytes(z["meta"])
    out = []
    for i, m in enumerate(raw.split(b"\n") if raw else []):
        a, b = z["offsets"][i], z["offsets"][i + 1]
        r = json.loads(m)
        r["d"] = [round(float(v), 5) for v in z["d"][a:b]]
        r["I"] = [round(float(v), 2) for v in z["I"][a:b]]
        r["hkl"] = z["hkl"][a:b, : z["hkl_len"][i]].tolist()
        out.append(r)
    return out

def _value(text: str, tag: str) -> str:
    """A CIF tag's value, single-line or on the next line(s) (quoted or ;-delimited)."""
    m = re.search(rf"^\s*{tag}\s*(.*)$", text, re.M)
    if not m:
        return ""
    v = m.group(1).strip()
    if not v:
        rest = text[m.end():].lstrip("\n")
        if rest.startswith(";"):
            v = rest[1:].split("\n;", 1)[0]
        else:
            v = rest.split("\n", 1)[0]
    v = " ".join(v.split()).strip("'\"")
    return "" if v in ("?", ".") else v


def wanted(path: Path) -> tuple[bool, str, str]:
    """(keep, reason or name, formula) from the CIF text alone, before any parsing."""
    if not path.stem.isdigit():
        return False, "not a COD number", ""
    if int(path.stem) >= AMCSD_FROM:
        return False, "AMCSD-derived", ""
    text = path.read_text(errors="replace")
    formula = _value(text, "_chemical_formula_sum")
    els = set(re.findall(r"[A-Z][a-z]?", formula))
    if {"C", "H"} <= els:
        return False, "organic", formula
    name = (_value(text, "_chemical_name_mineral") or _value(text, "_chemical_name_common")
            or _value(text, "_chemical_name_systematic"))
    return True, name, formula


def _alarm(*_):
    raise TimeoutError


def one(path: str):
    p = Path(path)
    try:
        keep, name, _ = wanted(p)
        if not keep:
            return None, name
        signal.signal(signal.SIGALRM, _alarm)
        signal.alarm(TIMEOUT)
        try:
            from pymatgen.io.cif import CifParser

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                s = CifParser(str(p)).parse_structures(primitive=False)[0]
            if len(s) > MAX_SITES:
                return None, "too many sites"
            r = reference_from_structure(s, ref_id=f"COD {p.stem}", source="Crystallography Open Database",
                                             licence="COD (open data, free for any use)")
        finally:
            signal.alarm(0)
        if name and len(name) < 80:
            r["name"] = f"{name} ({r['formula']})" if name.lower() != r["formula"].lower() else r["formula"]
        r["url"] = f"https://www.crystallography.net/cod/{p.stem}.html"
        keep_lines = [i for i, v in enumerate(r["I"]) if v >= MIN_I]
        for k in ("d", "I", "hkl"):
            r[k] = [r[k][i] for i in keep_lines]
        return r, ""
    except TimeoutError:
        return None, "timed out"
    except Exception as e:  # malformed, unparseable or missing-coordinate CIFs are common in a bulk download
        return None, f"error: {type(e).__name__}"


def shard(cif_dir: str, out: str, workers: int):
    files = sorted(str(p) for p in Path(cif_dir).rglob("*.cif"))
    refs, why = [], {}
    with Pool(workers, maxtasksperchild=200) as pool:
        for i, (r, reason) in enumerate(pool.imap_unordered(one, files, chunksize=8), 1):
            if r:
                refs.append(r)
            else:
                why[reason] = why.get(reason, 0) + 1
            if i % 2000 == 0:
                print(f"{i}/{len(files)} CIFs, {len(refs)} phases", flush=True)
    refs.sort(key=lambda r: int(r["id"].split()[-1]))
    save(refs, out)
    print(f"{len(files)} CIFs: {len(refs)} phases written to {out}; left out: {why}")


def merge(out: str, shards: list[str]):
    seen, refs = set(), []
    for p in shards:
        for r in load(p):
            if r["id"] not in seen:
                seen.add(r["id"])
                refs.append(r)
    refs.sort(key=lambda r: int(r["id"].split()[-1]))
    save(refs, out)
    print(f"{len(refs)} phases, {sum(len(r['d']) for r in refs)} lines written to {out}")


def plan(listing: str, max_slices: int = 240, min_per: int = 2500):
    counts = {}
    for line in open(listing, errors="replace"):
        path = line.rsplit(None, 1)[-1] if line.strip() else ""
        m = re.fullmatch(r"(\d)/(\d\d)/(\d\d)/(\d+)\.cif", path)
        if m and int(m.group(4)) < AMCSD_FROM:
            unit = "/".join(m.groups()[:3])
            counts[unit] = counts.get(unit, 0) + 1
    total = sum(counts.values())
    if not total:
        raise SystemExit("no COD CIF paths found in the listing")
    per = max(min_per, -(-total // max_slices))
    slices, cur, n = [], [], 0
    for unit in sorted(counts):
        cur.append(unit)
        n += counts[unit]
        if n >= per:
            slices.append(",".join(cur))
            cur, n = [], 0
    if cur:
        slices.append(",".join(cur))
    print(json.dumps(slices))
    import sys
    print(f"{total} CIFs below {AMCSD_FROM} in {len(counts)} folders -> {len(slices)} slices of about {per}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("shard")
    a.add_argument("cif_dir")
    a.add_argument("out")
    a.add_argument("--workers", type=int, default=4)
    pl = sub.add_parser("plan")
    pl.add_argument("listing")
    m = sub.add_parser("merge")
    m.add_argument("out")
    m.add_argument("shards", nargs="+")
    args = ap.parse_args()
    if args.cmd == "shard":
        shard(args.cif_dir, args.out, args.workers)
    elif args.cmd == "plan":
        plan(args.listing)
    else:
        merge(args.out, args.shards)


if __name__ == "__main__":
    main()
