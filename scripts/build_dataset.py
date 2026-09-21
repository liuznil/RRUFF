"""
Build the merged structure--chemistry--Raman dataset.

Synchronised with the revised ``structure_features.py`` (revision 2).

Design
------
1. AMCSD and RRUFF are matched by normalised mineral name.
2. One representative AMCSD structure and one representative Raman spectrum are chosen per matched
   mineral name.
3. Every field of ``StructureFeatures`` / ``RamanFeatures`` is propagated automatically, so new
   features never need a manual column list.
4. NAME COLLISIONS ARE RESOLVED EXPLICITLY.  Both dataclasses define ``ok``, ``error`` and
   ``mineral``.  The previous build did ``row.update(raman_row)``, which silently overwrote the
   structure-side ``ok`` / ``error`` / ``mineral`` (all structure diagnostics such as
   ``sg_mismatch_after_primitive`` or ``bond_failed`` were lost).  Here the ambiguous columns are
   renamed (``structure_ok``, ``structure_error``, ``raman_ok``, ``raman_error`` ...) and a check
   raises if any collision remains.
5. Every dropped mineral is logged with stage and reason (``dataset_failures.csv``).
6. Structure extraction (the slow part: symmetry + SOAP) can run in parallel: ``-j N``.

Usage
-----
    python build_dataset.py                 # project root = parent of this file's folder
    python build_dataset.py --root /path/to/project -j 8
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, fields, is_dataclass
from functools import lru_cache
from pathlib import Path

import pandas as pd

from raman_features import RamanFeatures, compute_raman_features
from structure_features import StructureFeatures, compute_structure_features

# -----------------------------------------------------------------------------
# Layout (expected)
#   <root>/data/AMCSD/cif/<Mineral>__<AMCSD_ID>.cif
#   <root>/data/raman/<quality>_<orientation>/<Mineral>__...Raman_Data_Processed...txt
#   <root>/results/
# -----------------------------------------------------------------------------
DEFAULT_ROOT = Path(__file__).resolve().parent.parent

# Higher-quality spectra first; within a tier unoriented before oriented (oriented single-crystal
# spectra are affected by polarisation selection rules).
RAMAN_TIER_ORDER = [
    ("excellent", "unoriented"),
    ("excellent", "oriented"),
    ("fair", "unoriented"),
    ("fair", "oriented"),
    ("unrated", "unoriented"),
    ("unrated", "oriented"),
    ("poor", "unoriented"),
    ("poor", "oriented"),
]

# Word boundaries: the old pattern "gpa|kbar" also matched inside unrelated words.
HIGH_PT_PATTERN = re.compile(
    r"\bgpa\b|\bkbar\b|high[ -]pressure|high[ -]temperature", re.IGNORECASE
)

# Columns that exist in BOTH feature classes (or are ambiguous in the merged table) are renamed.
STRUCT_RENAME = {"ok": "structure_ok", "error": "structure_error"}
RAMAN_RENAME = {
    "mineral": "rruff_mineral",
    "file_path": "raman_path",
    "quality": "raman_quality",
    "orientation": "raman_orientation",
    "wavelength": "raman_wavelength",
    "n_points": "raman_n_points",
    "snr": "raman_snr",
    "ok": "raman_ok",
    "error": "raman_error",
}

STRUCT_STAGE_CATEGORIES = ("parse_failed", "symmetry_failed", "group_theory_failed")


# -----------------------------------------------------------------------------
# Name normalisation and indexing
# -----------------------------------------------------------------------------
def norm(name: str) -> str:
    """Normalise a mineral name for AMCSD--RRUFF matching."""
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


def index_cifs(cif_dir: Path) -> dict[str, list[Path]]:
    """Index AMCSD CIFs by normalised mineral name (file name: <Mineral>__<AMCSD_ID>.cif)."""
    if not cif_dir.exists():
        raise FileNotFoundError(f"AMCSD CIF directory not found: {cif_dir}")
    by_name: dict[str, list[Path]] = defaultdict(list)
    for p in sorted(cif_dir.glob("*.cif")):
        by_name[norm(p.stem.rsplit("__", 1)[0])].append(p)
    return by_name


def index_raman(raman_root: Path) -> dict[str, dict[tuple[str, str], list[Path]]]:
    """Index processed RRUFF Raman spectra by normalised mineral name and (quality, orientation)."""
    by_name: dict[str, dict[tuple[str, str], list[Path]]] = defaultdict(lambda: defaultdict(list))
    n_dirs = 0
    for quality, orientation in RAMAN_TIER_ORDER:
        d = raman_root / f"{quality}_{orientation}"
        if not d.exists():
            continue
        n_dirs += 1
        for p in sorted(d.iterdir()):
            if p.is_file() and p.suffix.lower() == ".txt" and "Raman_Data_Processed" in p.name:
                by_name[norm(p.stem.split("__", 1)[0])][(quality, orientation)].append(p)
    if n_dirs == 0:
        raise FileNotFoundError(f"no <quality>_<orientation> folders found under {raman_root}")
    return by_name


# -----------------------------------------------------------------------------
# Representative selection
# -----------------------------------------------------------------------------
@lru_cache(maxsize=None)
def _is_high_pt(path_str: str) -> bool:
    try:
        text = Path(path_str).read_text(encoding="latin-1", errors="ignore")
    except OSError:
        return False
    return HIGH_PT_PATTERN.search(text) is not None


def _cif_priority_key(path: Path) -> tuple[int, int, str]:
    """Ambient-condition refinements first, then lower AMCSD id (deterministic)."""
    m = re.search(r"__(\d+)\.cif$", path.name, flags=re.IGNORECASE)
    return (1 if _is_high_pt(str(path)) else 0, int(m.group(1)) if m else 10**9, path.name)


def _structure_category(error: str) -> str:
    for cat in STRUCT_STAGE_CATEGORIES:
        if cat in error:
            return cat
    return "no_atoms" if not error else "other"


def pick_structure(paths: list[Path]):
    """Return (StructureFeatures | None, failure_category, failure_detail)."""
    first_err = first_cat = ""
    for p in sorted(paths, key=_cif_priority_key):
        feat = compute_structure_features(p)
        if feat.ok and feat.n_atom_primitive > 0:
            return feat, "", ""
        if not first_err:                 # report the failure of the highest-priority candidate
            first_err, first_cat = f"{p.name}: {feat.error[:200]}", _structure_category(feat.error)
    return None, (first_cat or "no_cif"), first_err


def _raman_priority_key(path: Path) -> tuple[int, str]:
    """Prefer 532 nm excitation, then deterministic file-name order."""
    return (0 if "__532__" in path.name else 1), path.name


def pick_raman(by_tier):
    """Return (RamanFeatures | None, failure_category, failure_detail)."""
    first_err = ""
    for tier in RAMAN_TIER_ORDER:
        for p in sorted(by_tier.get(tier, []), key=_raman_priority_key):
            feat = compute_raman_features(p, *tier)
            if feat.ok:
                return feat, "", ""
            first_err = first_err or f"{p.name}: {feat.error[:200]}"
    cat = re.match(r"[a-z_]+", first_err.split(": ", 1)[-1]) if first_err else None
    return None, (cat.group(0) if cat else "no_spectrum"), first_err


# -----------------------------------------------------------------------------
# Row construction
# -----------------------------------------------------------------------------
def _as_dict(obj) -> dict:
    if hasattr(obj, "to_dict"):
        return obj.to_dict()
    if is_dataclass(obj):
        return asdict(obj)
    raise TypeError(f"{type(obj).__name__} has neither to_dict() nor dataclass structure")


def merge_row(sfeat, rfeat, mineral_norm: str) -> dict:
    """Merge structure and Raman features; ambiguous names are renamed, never overwritten."""
    srow = {STRUCT_RENAME.get(k, k): v for k, v in _as_dict(sfeat).items()}
    rrow = {RAMAN_RENAME.get(k, k): v for k, v in _as_dict(rfeat).items()}
    clash = set(srow) & set(rrow)
    if clash:
        raise KeyError(f"unresolved column-name collision between structure and Raman: {sorted(clash)}")
    return {"mineral_norm": mineral_norm, **srow, **rrow}


def _process_one(task):
    """Worker: one mineral. Returns (key, row | None, (stage, category, detail) | None)."""
    key, cif_paths, raman_tiers = task
    sfeat, cat, detail = pick_structure(cif_paths)
    if sfeat is None:
        return key, None, ("structure", cat, detail)
    rfeat, cat, detail = pick_raman(raman_tiers)
    if rfeat is None:
        return key, None, ("raman", cat, detail)
    return key, merge_row(sfeat, rfeat, key), None


PREFERRED_COLUMNS = [
    # identity / matching
    "mineral", "mineral_norm", "rruff_mineral", "formula", "family", "amcsd_id", "rruff_id",
    "cif_path", "raman_path",
    # parsing / disorder / status
    "structure_ok", "structure_error", "cif_repaired", "n_unknown_labels",
    "is_ordered", "n_disordered_sites", "mixing_entropy",
    # symmetry
    "space_group_symbol", "space_group_number", "crystal_system", "point_group",
    "is_centrosymmetric", "symmetry_consistent", "n_sites_conventional", "n_atom_primitive",
    "z_primitive", "n_wyckoff_orbits",
    # group theory
    "n_vibrational_modes", "n_vib_distinct", "n_raman_active_modes", "n_ir_active_modes",
    "n_raman_ir_overlap", "n_silent_modes", "n_raman_dof", "n_ir_dof", "raman_mode_fraction",
    "degeneracy_fraction", "n_raman_irreps", "max_raman_irrep_mult", "raman_channel_entropy",
    "group_theory_status",
    # composition / local structure / SOAP
    "sigma_chi", "mass_contrast",
    "n_bonds", "bonds_per_atom", "mean_bond", "cv_bond", "bond_distortion", "mean_cn", "cn_std",
    "mean_delta_chi_bond", "mean_reduced_mass", "n_anion_group_types", "anion_group_o_fraction",
    "polymerization_index",
    "soap_intra_species_dispersion", "soap_n_env", "soap_min_env_sep",
    # Raman metadata / descriptors
    "raman_ok", "raman_error", "raman_quality", "raman_orientation", "raman_wavelength",
    "raman_n_points", "raman_snr",
    "n_peak", "s_peak", "gamma", "domega_median", "r_overlap",
    "n_low", "s_low", "gamma_low", "n_mid", "s_mid", "gamma_mid", "n_high", "s_high", "gamma_high",
    "w1_distance",
]


def _ordered_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Stable column order; columns not in PREFERRED_COLUMNS are kept and appended."""
    head = [c for c in PREFERRED_COLUMNS if c in df.columns]
    return df[head + [c for c in df.columns if c not in head]]


def expected_columns() -> set[str]:
    """Columns the merged table must contain, derived from the dataclasses (no hand-kept list)."""
    s = {STRUCT_RENAME.get(f.name, f.name) for f in fields(StructureFeatures)}
    r = {RAMAN_RENAME.get(f.name, f.name) for f in fields(RamanFeatures)}
    return s | r | {"mineral_norm"}


# -----------------------------------------------------------------------------
# Main build
# -----------------------------------------------------------------------------
def build(root: Path = DEFAULT_ROOT, n_jobs: int = 1) -> pd.DataFrame:
    cif_dir = root / "data" / "AMCSD" / "cif"
    raman_root = root / "data" / "raman"
    out_dir = root / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path, fail_path = out_dir / "dataset.csv", out_dir / "dataset_failures.csv"

    cif_index, raman_index = index_cifs(cif_dir), index_raman(raman_root)
    common = sorted(set(cif_index) & set(raman_index))
    print(f"AMCSD mineral names       : {len(cif_index)}")
    print(f"RRUFF Raman mineral names : {len(raman_index)}")
    print(f"Matched mineral names     : {len(common)}")

    tasks = [(k, cif_index[k], {t: list(v) for t, v in raman_index[k].items()}) for k in common]
    if n_jobs > 1:
        with ProcessPoolExecutor(max_workers=n_jobs) as ex:
            results = ex.map(_process_one, tasks, chunksize=4)      # order-preserving
            done = _collect(results, len(tasks))
    else:
        done = _collect(map(_process_one, tasks), len(tasks))
    rows, failures = done

    df = pd.DataFrame(rows)
    if not df.empty:
        df = _ordered_columns(df)
    df.to_csv(out_path, index=False)
    fail_df = pd.DataFrame(failures, columns=["mineral_norm", "stage", "category", "detail"])
    fail_df.to_csv(fail_path, index=False)

    _report(df, fail_df, len(common), out_path, fail_path)
    return df


def _collect(results, total: int):
    rows, failures = [], []
    for i, (key, row, fail) in enumerate(results, start=1):
        if row is not None:
            rows.append(row)
        else:
            failures.append((key, *fail))
        if i % 100 == 0 or i == total:
            print(f"  ... {i}/{total} processed, {len(rows)} kept, {len(failures)} dropped", flush=True)
    return rows, failures


def _report(df: pd.DataFrame, fail_df: pd.DataFrame, n_matched: int, out_path: Path, fail_path: Path):
    print("\n" + "=" * 72 + "\nDataset build completed\n" + "=" * 72)
    print(f"Matched mineral names : {n_matched}")
    print(f"Rows kept             : {len(df)}")
    print(f"Rows dropped          : {len(fail_df)}")
    if len(fail_df):
        for (stage, cat), n in Counter(zip(fail_df["stage"], fail_df["category"])).most_common():
            print(f"    {stage:9s} {cat:22s} {n}")
    print(f"Output                : {out_path}")
    print(f"Failure log           : {fail_path}")

    if df.empty:
        print("WARNING: dataset is empty.")
        return
    missing = sorted(expected_columns() - set(df.columns))
    print("WARNING: missing columns: " + ", ".join(missing) if missing
          else "All structure and Raman feature columns are present.")
    if "symmetry_consistent" in df:
        n_bad = int((~df["symmetry_consistent"].astype(bool)).sum())
        if n_bad:
            print(f"NOTE: {n_bad} rows have symmetry_consistent=False "
                  "(group-theory features describe a different space group than space_group_*).")
    if "structure_error" in df:
        n_err = int(df["structure_error"].fillna("").astype(str).str.len().gt(0).sum())
        print(f"NOTE: {n_err} kept rows carry a non-empty structure_error (partial descriptors).")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Build the merged structure-Raman dataset")
    ap.add_argument("--root", type=Path, default=DEFAULT_ROOT, help="project root (contains data/)")
    ap.add_argument("-j", "--jobs", type=int, default=1, help="parallel worker processes")
    a = ap.parse_args(argv)
    build(a.root, a.jobs)
    return 0


if __name__ == "__main__":
    sys.exit(main())