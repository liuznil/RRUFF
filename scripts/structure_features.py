"""
Structural / symmetry feature extraction from AMCSD CIF files.

Research logic implemented here
-------------------------------
  crystal structure  -> "latent mode space"          (N atoms  -> 3N-3 optical dof)
  symmetry / motifs  -> how modes split & cluster    (irreps, degeneracy, multiplicities,
                                                      inequivalent sites, anion groups)
  Raman experiment   -> finite number of peaks       (Raman-active *distinct* frequencies)
  chemistry / local  -> peak position, spacing and   (bond length / distortion / reduced mass,
  structure             resolvability                 electronegativity contrast, SOAP
                                                      environment dispersion)

Main changes w.r.t. the previous version (see the accompanying review for details)
---------------------------------------------------------------------------------
* Group theory: the hand-typed 32-point-group character tables (17 of 32 were invalid; the
  "(det, trace) identifies the class" assumption is false) were removed. Irreducible characters are
  now *computed* from the group multiplication table (Burnside-Dixon); Raman / IR activity is
  derived from Sym^2(V) / V. Verified against quartz, rutile, forsterite, NaCl, diamond.
* Fixed-atom test uses a Cartesian tolerance (A) instead of 1e-4 in fractional coordinates.
* Mode counts distinguish "distinct frequencies" (what a spectrum can show) from degrees of freedom.
* `cv_bond` is now a true coefficient of variation; the old quantity is `bond_distortion`.
* `soap_distortion` (variance of raw SOAP components, dominated by chemical contrast and by the
  number of species) replaced by symmetry-inequivalent-environment measures on normalised SOAP.
* Legacy-CIF label repair no longer maps Ba/Pb/Sr/Sn/Sb/Bi/Be/Se/Hg/... onto the wrong element.
* No global warning suppression, no print(), lazy heavy imports, pymatgen API-version shims.

Revision 2 (verification pass)
------------------------------
* Bonds: a bond i->j now requires j to be "anion-like" (chi within ANION_CHI_WINDOW of the most
  electronegative species present). The old rule "chi_j > chi_i" also counted cation-cation contacts
  (Mg->Si, Ba->S, Ca->Al ...), inflating mean_cn / cv_bond / bond_distortion / mean_delta_chi_bond
  (forsterite: mean_cn 7.67 -> 5.33).
* Z and formula are computed from the site-based (majority-species) cell, so disordered and
  non-stoichiometric structures no longer give Z<1 or formulas such as "Fe3.68O4".
* New field `symmetry_consistent`: False when the primitive cell's space group differs from the one
  found for the CIF cell (group-theory features then describe a different group than `space_group_*`).
* Legacy-CIF repair ignores '#' comment lines (they used to be counted as unknown atom labels).

Only numpy is needed to import this module; pymatgen (and dscribe for SOAP) are imported lazily.
"""
from __future__ import annotations

import argparse
import csv
import logging
import math
import re
import warnings
from dataclasses import asdict, dataclass, fields
from functools import lru_cache
from pathlib import Path

import numpy as np

log = logging.getLogger(__name__)

# --------------------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------------------
SYMPREC = 0.1                 # A, spglib tolerance (experimental CIFs)
ANGLE_TOLERANCE = 5.0         # deg
FIXED_ATOM_TOL = 2.0 * SYMPREC  # A, atom-maps-onto-itself tolerance (Cartesian!)

BOND_TOLERANCE = 1.3          # bond if d <= 1.3 * (r_cov_i + r_cov_j)
MAX_SEARCH_RADIUS = 6.0       # A
GAP_RATIO_THRESHOLD = 1.3     # first coordination shell ends at first gap >= 30 %
ANION_CHI_WINDOW = 0.5        # a species is "anion-like" if chi >= chi_max(structure) - window
DEFAULT_RADIUS = 1.5          # A, used only if no radius is available at all

SOAP_R_CUT = 5.0
SOAP_N_MAX = 6
SOAP_L_MAX = 4
SOAP_ENV_TOL = 0.01           # distance between unit-normalised SOAP vectors: "same environment"
MAX_ATOMS_FOR_SOAP = 500      # primitive-cell atoms

# Oxyanion / network formers: a (former, n_O) pair with n_O in {3, 4} is an "anion group"
ANION_FORMERS = frozenset({"B", "C", "N", "Si", "P", "S", "V", "Cr", "As", "Se", "Mo", "W", "Ge"})


# --------------------------------------------------------------------------------------
# Element data (lazy, cached)
# --------------------------------------------------------------------------------------
@lru_cache(maxsize=None)
def _element_data(symbol: str) -> tuple[float, float, float]:
    """(Pauling electronegativity, atomic mass, covalent radius in A) with safe fallbacks."""
    from pymatgen.core.periodic_table import Element

    chi, mass, rad = 2.0, 16.0, None
    try:
        el = Element(symbol)
        x = el.X
        if x is not None and not (isinstance(x, float) and math.isnan(x)):
            chi = float(x)
        mass = float(el.atomic_mass)
        rad = el.atomic_radius        # empirical radius, only a fallback
        rad = float(rad) if rad is not None else None
    except Exception:                 # dummy species "X", etc.
        pass
    try:
        from pymatgen.analysis.molecule_structure_comparator import CovalentRadius
        r_cov = CovalentRadius.radius.get(symbol)
        if r_cov is not None:
            rad = float(r_cov)
    except Exception:
        pass
    return chi, mass, (rad if rad else DEFAULT_RADIUS)


def _site_symbols(site):
    """[(symbol, occupancy), ...] for a (possibly disordered) site."""
    return [(sp.symbol, float(occ)) for sp, occ in site.species.items()]


def _site_average(site, idx: int) -> float:
    """Occupancy-weighted average of element property idx (0: chi, 1: mass, 2: r_cov)."""
    tot = val = 0.0
    for sym, occ in _site_symbols(site):
        val += occ * _element_data(sym)[idx]
        tot += occ
    return val / tot if tot > 0 else _element_data("X")[idx]


# --------------------------------------------------------------------------------------
# CIF loading
# --------------------------------------------------------------------------------------
_ELEMENTS = frozenset(
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se "
    "Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb "
    "Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm".split()
)


def _label_to_element(label: str) -> str | None:
    """'Ba1' -> 'Ba', 'OH2' -> 'O', 'Wat1' -> 'O', 'Sb' -> 'Sb'. None if unrecognised."""
    lab = label.strip("'\"")
    if lab.lower().startswith("wat"):
        return "O"
    if lab.upper() == "D":
        return "H"
    m = re.match(r"[A-Za-z]+", lab)
    if not m:
        return None
    alpha = m.group(0)
    if len(alpha) >= 2 and alpha[1].islower():
        two = alpha[0].upper() + alpha[1]
        if two in _ELEMENTS:
            return two
    one = alpha[0].upper()
    return one if one in _ELEMENTS else None


def _repair_legacy_cif(text: str) -> tuple[str, int]:
    """Repair old AMCSD-style CIFs. Returns (text, number of atom labels that could not be resolved).

    The previous version used a hand-written alternation regex in which e.g. "Se" was shadowed by "S"
    and Ba, Pb, Sr, Sn, Sb, Bi, Be, Hg, Hf, Ho ... silently became B, P, S, S, S, B, B, H, H, H.
    """
    repaired, in_atom_loop, n_unknown = [], False, 0
    for raw_line in text.splitlines():
        stripped = raw_line.strip()
        if stripped.lower() == "_space_group_symop_operation_xyz":
            raw_line = raw_line.replace("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz")
            stripped = raw_line.strip()
        if stripped.upper() == "END":
            continue
        if stripped.lower() == "_atom_site_label":
            in_atom_loop = True
        elif in_atom_loop and stripped.lower() == "loop_":
            in_atom_loop = False
        elif in_atom_loop and stripped and not stripped.startswith(("_", "#")):
            fields_ = stripped.split()
            if len(fields_) >= 4:
                el = _label_to_element(fields_[0])
                if el is None:
                    n_unknown += 1
                    el = "X"
                fields_[0] = el
                indent = raw_line[: len(raw_line) - len(raw_line.lstrip())]
                raw_line = indent + " ".join(fields_)
        repaired.append(raw_line)
    return "\n".join(repaired) + "\n", n_unknown


def _cif_structures(parser):
    """pymatgen renamed get_structures -> parse_structures (2024); support both."""
    fn = getattr(parser, "parse_structures", None) or parser.get_structures
    return fn(primitive=False)


def load_structure(cif_path: Path):
    """Returns (structure, flags). flags: repaired (bool), n_unknown_labels (int). Raises on failure."""
    from pymatgen.core import Structure
    from pymatgen.io.cif import CifParser

    flags = {"repaired": False, "n_unknown_labels": 0}
    errors = []
    try:
        return Structure.from_file(str(cif_path)), flags
    except Exception as e:
        errors.append(f"parse: {e}")
    try:
        structures = _cif_structures(CifParser(str(cif_path), occupancy_tolerance=1.0, check_cif=False))
        if structures:
            return structures[0], flags
        errors.append("fallback: no structures")
    except Exception as e:
        errors.append(f"fallback: {e}")
    try:
        text, n_unknown = _repair_legacy_cif(cif_path.read_text(encoding="latin-1", errors="replace"))
        from_str = getattr(CifParser, "from_str", None) or CifParser.from_string
        structures = _cif_structures(from_str(text, occupancy_tolerance=1.0, check_cif=False))
        if structures:
            flags.update(repaired=True, n_unknown_labels=n_unknown)
            return structures[0], flags
        errors.append("repair: no structures")
    except Exception as e:
        errors.append(f"repair: {e}")
    raise ValueError(" | ".join(errors))


# --------------------------------------------------------------------------------------
# Disorder
# --------------------------------------------------------------------------------------
def _make_ordered_structure(structure):
    """Majority-species ordered copy (deterministic tie-break by symbol), occupancy forced to 1.
    Only used for SOAP; symmetry and bond analysis use the original (average) structure."""
    from pymatgen.core import Structure

    species = []
    for site in structure:
        items = sorted(site.species.items(), key=lambda kv: (-kv[1], kv[0].symbol))
        species.append(items[0][0])
    return Structure(structure.lattice, species, structure.frac_coords, coords_are_cartesian=False)


def formula_units(structure) -> float:
    """Number of formula units in `structure`, from the site-based (majority-species) composition.
    Occupancy-weighted compositions have non-integer amounts, for which pymatgen's reduced formula
    is meaningless (Fe0.92O -> "Fe3.68O4"); the majority-species cell always has integer counts."""
    ref = structure if structure.is_ordered else _make_ordered_structure(structure)
    return float(ref.composition.get_reduced_composition_and_factor()[1])


def disorder_features(structure) -> dict:
    n_dis, ent = 0, []
    for site in structure:
        occ = [float(o) for o in site.species.values()]
        tot = sum(occ)
        if len(occ) > 1 or tot < 1 - 1e-6:
            n_dis += 1
        p = [o for o in occ if o > 0] + ([1.0 - tot] if tot < 1 - 1e-6 else [])
        ent.append(-sum(o * math.log(o) for o in p))
    return {"is_ordered": n_dis == 0, "n_disordered_sites": n_dis,
            "mixing_entropy": float(np.mean(ent)) if ent else 0.0}


# --------------------------------------------------------------------------------------
# Composition family (coarse, composition-only; see anion-group features for structure-based info)
# --------------------------------------------------------------------------------------
_HALOGENS = {"F", "Cl", "Br", "I"}
_CHALCOGENS_NO_O = {"S", "Se", "Te"}


def classify_family(comp) -> str:
    """comp: pymatgen Composition or {symbol: amount}. Coarse heuristic - use as a categorical
    grouping only. Thresholds are O:X ratios needed for X to be the centre of an oxyanion."""
    amt = comp.get_el_amt_dict() if hasattr(comp, "get_el_amt_dict") else dict(comp)
    els = {k for k, v in amt.items() if v > 0}
    n = lambda s: float(amt.get(s, 0.0))

    if els <= {"O", "H"}:
        return "oxide_hydroxide" if "O" in els else "native_element"
    if len(els - {"H"}) == 1 and "O" not in els:
        return "native_element"
    has_O = n("O") > 0
    if not has_O:
        if els & _CHALCOGENS_NO_O:
            return "sulfide_selenide_telluride"
        if els & _HALOGENS:
            return "halide"
        return "other"

    nO = n("O")
    ratio = lambda *xs: nO / max(sum(n(x) for x in xs), 1e-9)
    if n("Si") > 0:
        return "silicate"
    if n("C") > 0 and n("P") == 0 and n("S") <= n("C") and ratio("C") >= 3.0 - 1e-9:
        return "carbonate"
    if n("S") > 0 and ratio("S") >= 3.5 - 1e-9:
        return "sulfate"
    if n("Cr") > 0 and ratio("Cr") >= 3.5 - 1e-9:
        return "chromate"
    if (n("P") + n("As") + n("V")) > 0 and ratio("P", "As", "V") >= 3.0 - 1e-9:
        return "phosphate_arsenate_vanadate"
    if n("B") > 0:
        return "borate"
    if n("N") > 0 and n("C") == 0 and ratio("N") >= 3.0 - 1e-9:
        return "nitrate"
    return "oxide_hydroxide"


# --------------------------------------------------------------------------------------
# Group theory (table-free)
# --------------------------------------------------------------------------------------
def _ikey(m) -> tuple:
    return tuple(int(x) for x in np.asarray(m).ravel())


def _group_tables(rots):
    n = len(rots)
    idx = {_ikey(r): i for i, r in enumerate(rots)}
    if len(idx) != n:
        raise ValueError("duplicate rotation matrices (cell is not primitive)")
    mul = np.empty((n, n), dtype=int)
    for a in range(n):
        for b in range(n):
            k = _ikey(rots[a] @ rots[b])
            if k not in idx:
                raise ValueError("rotation set is not closed under multiplication")
            mul[a, b] = idx[k]
    e = idx[_ikey(np.eye(3, dtype=int))]
    inv = np.array([int(np.flatnonzero(mul[a] == e)[0]) for a in range(n)])
    return mul, inv, e


def _conjugacy_classes(mul, inv):
    n = len(mul)
    cls_of = -np.ones(n, dtype=int)
    classes = []
    for g in range(n):
        if cls_of[g] >= 0:
            continue
        members = sorted({int(mul[mul[h, g], inv[h]]) for h in range(n)})
        for m in members:
            cls_of[m] = len(classes)
        classes.append(members)
    return cls_of, classes


def irreducible_characters(rots) -> dict:
    """Irreducible characters of a finite group given as integer matrices (Burnside-Dixon).

    Complex characters are returned as-is; pairing of complex-conjugate irreps is done by the caller.
    """
    rots = [np.asarray(r, dtype=int) for r in rots]
    n = len(rots)
    mul, inv, e = _group_tables(rots)
    cls_of, classes = _conjugacy_classes(mul, inv)
    nc = len(classes)
    sizes = np.array([len(c) for c in classes], dtype=float)
    # class-multiplication constants a_ijk = #{(x,y) in C_i x C_j : xy = z_k}
    A = np.zeros((nc, nc, nc))
    for k in range(nc):
        z = classes[k][0]
        for i in range(nc):
            for x in classes[i]:
                A[i, cls_of[mul[inv[x], z]], k] += 1
    ic = int(cls_of[e])
    rng = np.random.default_rng(12345)
    for _ in range(6):                     # a random combination has non-degenerate eigenvalues
        _, v = np.linalg.eig(np.tensordot(rng.normal(size=nc), A, axes=1))
        try:
            chars = []
            for m in range(nc):
                om = v[:, m] / v[ic, m]
                d = math.sqrt(n / np.sum(np.abs(om) ** 2 / sizes))
                if abs(d - round(d)) > 1e-6:
                    raise ValueError("non-integer irrep dimension")
                chars.append(om * round(d) / sizes)
            chars = np.array(chars)
            gram = (np.conj(chars) * sizes) @ chars.T / n
            if np.allclose(gram, np.eye(nc), atol=1e-6):
                return dict(n=n, mul=mul, inv=inv, e=e, cls_of=cls_of, classes=classes,
                            sizes=sizes, chars=chars)
        except ValueError:
            continue
    raise ValueError("could not obtain orthonormal irreducible characters")


@dataclass
class GroupTheoryFeatures:
    point_group: str = ""
    n_vibrational_modes: int = 0        # 3N-3, degrees of freedom of the primitive cell
    n_vib_distinct: int = 0             # distinct frequencies (degenerate irreps counted once)
    n_raman_active_modes: int = 0       # distinct Raman-active frequencies (what a spectrum can show)
    n_ir_active_modes: int = 0
    n_raman_ir_overlap: int = 0
    n_silent_modes: int = 0
    n_raman_dof: int = 0                # counted with degeneracy
    n_ir_dof: int = 0
    raman_mode_fraction: float = float("nan")   # n_raman_active / n_vib_distinct
    degeneracy_fraction: float = float("nan")   # 1 - n_vib_distinct / (3N-3)
    n_raman_irreps: int = 0             # number of populated Raman symmetry channels
    max_raman_irrep_mult: int = 0       # most modes in one Raman irrep (same-symmetry modes repel)
    raman_channel_entropy: float = float("nan")
    is_centrosymmetric: bool = False
    status: str = "not_run"
    error: str = ""


def analyse_modes(rots, n_fixed) -> dict:
    """Gamma-point mode analysis of a PRIMITIVE cell.

    rots    : (n,3,3) integer rotation matrices (fractional basis) of the factor group
    n_fixed : (n,) number of atoms mapped onto themselves (mod lattice) by each operation

    chi_disp(g) = n_fixed(g) * tr(R);  Gamma_vib = Gamma_disp - Gamma_acoustic (= V, chi = tr R)
    Raman-active irreps occur in Sym^2(V), IR-active irreps occur in V.
    Complex-conjugate irrep pairs (time reversal) are merged into one doubly degenerate level.
    """
    rots = np.asarray(np.rint(rots), dtype=int)
    G = irreducible_characters(list(rots))
    n, mul, cls_of, sizes, chars = G["n"], G["mul"], G["cls_of"], G["sizes"], G["chars"]

    tr = np.array([np.trace(r) for r in rots], dtype=float)          # chi_V; similarity-invariant
    tr_sq = tr[[mul[g, g] for g in range(n)]]
    chi_disp = np.asarray(n_fixed, dtype=float) * tr
    chi_vib = chi_disp - tr
    chi_sym2 = 0.5 * (tr ** 2 + tr_sq)

    def to_class_function(f):
        out = np.zeros(len(sizes))
        for c, members in enumerate(G["classes"]):
            vals = f[members]
            if np.ptp(vals) > 1e-9:
                raise ValueError("character is not constant on a conjugacy class "
                                 "(atom mapping inconsistent with the symmetry operations)")
            out[c] = vals[0]
        return out

    def multiplicity(f):
        m = (np.conj(chars) * sizes * to_class_function(f)).sum(axis=1) / n
        if np.abs(m.imag).max() > 1e-6 or np.abs(m.real - np.rint(m.real)).max() > 1e-6:
            raise ValueError("non-integer irrep multiplicities")
        return np.rint(m.real).astype(int)

    m_vib, m_V, m_S2 = multiplicity(chi_vib), multiplicity(tr), multiplicity(chi_sym2)
    if (m_vib < 0).any():
        raise ValueError("negative irrep multiplicity")

    dims = np.rint(chars[:, cls_of[G["e"]]].real).astype(int)
    used = [False] * len(chars)
    phys = []
    for i in range(len(chars)):
        if used[i]:
            continue
        used[i] = True
        paired = False
        if np.abs(chars[i].imag).max() > 1e-8:
            for j in range(i + 1, len(chars)):
                if not used[j] and np.allclose(chars[j], np.conj(chars[i]), atol=1e-6):
                    used[j] = paired = True
                    break
            if not paired:
                raise ValueError("unpaired complex irrep")
        phys.append(dict(dim=int(dims[i]) * (2 if paired else 1), mult=int(m_vib[i]),
                         raman=bool(m_S2[i] > 0), ir=bool(m_V[i] > 0)))

    n3 = int(round(chi_disp[G["e"]]))
    total_dof = sum(p["mult"] * p["dim"] for p in phys)
    if total_dof != n3 - 3:
        raise ValueError(f"mode count mismatch: {total_dof} != 3N-3 = {n3 - 3}")

    ram = [p for p in phys if p["raman"] and p["mult"] > 0]
    w = np.array([p["mult"] for p in ram], dtype=float)
    entropy = float(-(w / w.sum() * np.log(w / w.sum())).sum()) if len(ram) else 0.0
    S = lambda cond, key="mult": sum(p["mult"] * (p["dim"] if key == "dof" else 1) for p in phys if cond(p))
    return dict(
        n_ops=n, n_vib_dof=total_dof, n_vib_distinct=S(lambda p: True),
        n_raman=S(lambda p: p["raman"]), n_ir=S(lambda p: p["ir"]),
        n_overlap=S(lambda p: p["raman"] and p["ir"]),
        n_silent=S(lambda p: not p["raman"] and not p["ir"]),
        n_raman_dof=S(lambda p: p["raman"], "dof"), n_ir_dof=S(lambda p: p["ir"], "dof"),
        n_raman_irreps=len(ram), max_raman_mult=int(w.max()) if len(ram) else 0,
        raman_entropy=entropy,
        centrosymmetric=any(_ikey(r) == _ikey(-np.eye(3, dtype=int)) for r in rots),
        irreps=phys,
    )


def count_fixed_atoms(rots, trans, frac_coords, lattice_matrix, tol: float = FIXED_ATOM_TOL):
    """Atoms mapped onto themselves (mod lattice) by each {R|t}. `tol` is a CARTESIAN distance (A).
    (The old code used 1e-4 in fractional units, which fails for any real, slightly noisy CIF.)"""
    frac = np.asarray(frac_coords, float)
    lat = np.asarray(lattice_matrix, float)
    out = np.zeros(len(rots), dtype=int)
    for k, (R, t) in enumerate(zip(rots, trans)):
        d = frac @ np.asarray(R, float).T + np.asarray(t, float) - frac
        d -= np.rint(d)
        out[k] = int((np.linalg.norm(d @ lat, axis=1) < tol).sum())
    return out


def compute_group_theory_features(structure, point_group: str, symmetry_operations,
                                  tol: float = FIXED_ATOM_TOL) -> GroupTheoryFeatures:
    """structure: PRIMITIVE pymatgen Structure; symmetry_operations: SpacegroupAnalyzer of that same
    structure, .get_symmetry_operations(cartesian=False). `point_group` is stored as metadata only."""
    res = GroupTheoryFeatures(point_group=point_group)
    res.n_vibrational_modes = max(0, 3 * len(structure) - 3)
    try:
        rots = np.array([op.rotation_matrix for op in symmetry_operations], dtype=float)
        if np.abs(rots - np.rint(rots)).max() > 1e-6:
            raise ValueError("non-integer fractional rotation matrices")
        rots = np.rint(rots).astype(int)
        trans = np.array([op.translation_vector for op in symmetry_operations], dtype=float)
        nfix = count_fixed_atoms(rots, trans, structure.frac_coords, structure.lattice.matrix, tol)
        m = analyse_modes(rots, nfix)
    except Exception as exc:
        res.status, res.error = "failed", f"{type(exc).__name__}: {exc}"
        return res
    res.n_vib_distinct = m["n_vib_distinct"]
    res.n_raman_active_modes, res.n_ir_active_modes = m["n_raman"], m["n_ir"]
    res.n_raman_ir_overlap, res.n_silent_modes = m["n_overlap"], m["n_silent"]
    res.n_raman_dof, res.n_ir_dof = m["n_raman_dof"], m["n_ir_dof"]
    if m["n_vib_distinct"] > 0:
        res.raman_mode_fraction = m["n_raman"] / m["n_vib_distinct"]
    if m["n_vib_dof"] > 0:
        res.degeneracy_fraction = 1.0 - m["n_vib_distinct"] / m["n_vib_dof"]
    res.n_raman_irreps, res.max_raman_irrep_mult = m["n_raman_irreps"], m["max_raman_mult"]
    res.raman_channel_entropy = m["raman_entropy"]
    res.is_centrosymmetric = m["centrosymmetric"]
    res.status = "ok"
    return res


# --------------------------------------------------------------------------------------
# Bonds / local coordination (pure function on neighbour lists -> testable without pymatgen)
# --------------------------------------------------------------------------------------
def _first_shell(dists):
    """Keep distances up to the first gap of >= GAP_RATIO_THRESHOLD (sorted)."""
    s = sorted(dists)
    for k in range(len(s) - 1):
        if s[k + 1] >= GAP_RATIO_THRESHOLD * s[k]:
            return s[: k + 1]
    return s


def analyse_bonds(neighbors, symbols, chis, masses, radii, homopolar: bool = False) -> dict:
    """neighbors[i] = [(j, distance), ...]; symbols[i] = dominant element symbol of site i.

    A bond i->j is counted when d <= BOND_TOLERANCE*(r_i+r_j), j is more electronegative than i AND j is
    anion-like (chi_j >= max(chi) - ANION_CHI_WINDOW).  The second condition removes cation-cation
    contacts such as Mg->Si or Ba->S that pass the plain electronegativity test.
    For single-element structures all contacts count.
    """
    nsite = len(neighbors)
    chi_max = max(chis) if len(chis) else 0.0
    anion_like = [c >= chi_max - ANION_CHI_WINDOW - 1e-9 for c in chis]
    site_bond_dist, site_cv, site_cn, all_d, dchi, mu = [], [], [], [], [], []
    o_former_count = {}                         # O site -> number of former neighbours
    group_types = set()
    for i, nbrs in enumerate(neighbors):
        cand = []
        for j, d in nbrs:
            ok = (chis[j] >= chis[i] - 1e-9) if homopolar else (anion_like[j] and chis[j] > chis[i] + 1e-9)
            if ok and d <= (radii[i] + radii[j]) * BOND_TOLERANCE:
                cand.append((d, j))
        cand = sorted(cand)
        keep = _first_shell([d for d, _ in cand])
        cand = cand[: len(keep)]
        # anion-group detection: former surrounded only by O, CN 3 or 4 (independent of the stats below)
        if symbols[i] in ANION_FORMERS and len(cand) in (3, 4) and all(symbols[j] == "O" for _, j in cand):
            group_types.add((symbols[i], len(cand)))
            for _, j in cand:
                o_former_count[j] = o_former_count.get(j, 0) + 1
        if len(cand) < 2:          # terminal contacts (e.g. O-H) are not coordination polyhedra
            continue
        ds = np.array([d for d, _ in cand])
        site_cn.append(len(ds))
        lm = ds.mean()
        site_bond_dist.append(float(np.mean(np.abs(ds - lm) / lm)))       # Baur-type distortion
        site_cv.append(float(ds.std() / lm))                              # true CV
        all_d.extend(ds.tolist())
        for d, j in cand:
            dchi.append(abs(chis[j] - chis[i]))
            mu.append(masses[i] * masses[j] / (masses[i] + masses[j]))
    out = dict(n_bonds=len(all_d), bonds_per_atom=len(all_d) / max(nsite, 1))
    nan = float("nan")
    out.update(mean_bond=float(np.mean(all_d)) if all_d else nan,
               cv_bond=float(np.mean(site_cv)) if site_cv else nan,
               bond_distortion=float(np.mean(site_bond_dist)) if site_bond_dist else nan,
               mean_cn=float(np.mean(site_cn)) if site_cn else nan,
               cn_std=float(np.std(site_cn)) if site_cn else nan,
               mean_delta_chi_bond=float(np.mean(dchi)) if dchi else nan,
               mean_reduced_mass=float(np.mean(mu)) if mu else nan)
    n_O = sum(1 for s in symbols if s == "O")
    out["n_anion_group_types"] = len(group_types)
    out["anion_group_o_fraction"] = (len(o_former_count) / n_O) if n_O else nan
    out["polymerization_index"] = (float(np.mean([c - 1 for c in o_former_count.values()]))
                                   if o_former_count else nan)   # 0: isolated groups, 1: all O bridge two
    return out


def bond_features(structure) -> dict:
    """pymatgen adapter for analyse_bonds (uses the ORIGINAL structure, i.e. experimental geometry)."""
    radii = [_site_average(s, 2) for s in structure]
    chis = [_site_average(s, 0) for s in structure]
    masses = [_site_average(s, 1) for s in structure]
    symbols = [max(_site_symbols(s), key=lambda t: t[1])[0] for s in structure]
    search_r = min(MAX_SEARCH_RADIUS, 2 * max(radii) * BOND_TOLERANCE + 0.5)
    all_nn = structure.get_all_neighbors(r=search_r, include_index=True)
    neighbors = [[(nb.index, nb.nn_distance) for nb in nbrs] for nbrs in all_nn]
    homopolar = len({s for s in symbols}) == 1
    return analyse_bonds(neighbors, symbols, chis, masses, radii, homopolar)


# --------------------------------------------------------------------------------------
# SOAP (local-environment heterogeneity, symmetry-aware)
# --------------------------------------------------------------------------------------
def soap_environment_stats(vectors, symbols, tol: float = SOAP_ENV_TOL) -> dict:
    """Statistics on unit-normalised SOAP vectors of all atoms of a cell.

    intra_species_dispersion : mean squared distance of an atom's SOAP vector to the centroid of the
        SOAP vectors of atoms of the SAME element (0 if all atoms of an element are equivalent).
        Unlike the variance of raw SOAP components it does not depend on the number of species /
        vector length, and it ignores the trivially large contrast between different elements.
    n_env : number of distinct local environments (same element, distance > tol).
    min_env_sep : smallest distance between two distinct environments of the same element
        (small => nearly degenerate, hard-to-resolve peaks). NaN if every element has one environment.
    """
    V = np.asarray(vectors, dtype=float)
    norms = np.linalg.norm(V, axis=1, keepdims=True)
    V = V / np.where(norms > 1e-12, norms, 1.0)
    total, n_env, seps = 0.0, 0, []
    for s in sorted(set(symbols)):
        sub = V[[i for i, x in enumerate(symbols) if x == s]]
        total += float(((sub - sub.mean(axis=0)) ** 2).sum())
        reps = []
        for v in sub:
            if not any(np.linalg.norm(v - r) <= tol for r in reps):
                reps.append(v)
        n_env += len(reps)
        for a in range(len(reps)):
            for b in range(a + 1, len(reps)):
                seps.append(float(np.linalg.norm(reps[a] - reps[b])))
    return dict(intra_species_dispersion=total / len(V), n_env=n_env,
                min_env_sep=min(seps) if seps else float("nan"))


def soap_features(prim) -> dict:
    """SOAP statistics on the (idealised) primitive cell; disorder -> majority species."""
    from dscribe.descriptors import SOAP
    from pymatgen.io.ase import AseAtomsAdaptor

    ordered = prim if prim.is_ordered else _make_ordered_structure(prim)
    if len(ordered) > MAX_ATOMS_FOR_SOAP:
        raise ValueError(f"{len(ordered)} atoms in primitive cell > {MAX_ATOMS_FOR_SOAP}")
    atoms = AseAtomsAdaptor.get_atoms(ordered)
    symbols = atoms.get_chemical_symbols()
    soap = SOAP(species=sorted(set(symbols)), periodic=True, r_cut=SOAP_R_CUT,
                n_max=SOAP_N_MAX, l_max=SOAP_L_MAX, average="off")
    try:
        vec = soap.create(atoms)
    except Exception:
        # Some dscribe builds may reject cells smaller than 2*r_cut: fall back to a supercell.
        # ASE keeps the first n atoms of repeat() identical to the original cell.
        reps = tuple(max(1, math.ceil((2 * SOAP_R_CUT + 0.5) / L)) for L in atoms.cell.lengths())
        big = atoms.repeat(reps)
        if len(big) > 4 * MAX_ATOMS_FOR_SOAP:
            raise
        vec = soap.create(big, centers=list(range(len(atoms))))
    if vec.ndim != 2 or vec.shape[0] != len(atoms) or not np.all(np.isfinite(vec)):
        raise ValueError(f"invalid SOAP output {getattr(vec, 'shape', None)}")
    return soap_environment_stats(vec, symbols)


# --------------------------------------------------------------------------------------
# Feature container
# --------------------------------------------------------------------------------------
_NAN = float("nan")


@dataclass
class StructureFeatures:
    # identity
    mineral: str
    amcsd_id: str
    cif_path: str
    formula: str
    family: str = "other"
    # parsing / disorder
    cif_repaired: bool = False
    n_unknown_labels: int = 0
    is_ordered: bool = True
    n_disordered_sites: int = 0
    mixing_entropy: float = _NAN
    # symmetry
    space_group_symbol: str = ""
    space_group_number: int = 0
    crystal_system: str = ""
    point_group: str = ""
    is_centrosymmetric: bool = False
    n_sites_conventional: int = 0
    n_atom_primitive: int = 0
    z_primitive: float = _NAN
    n_wyckoff_orbits: int = 0
    symmetry_consistent: bool = True    # False: primitive-cell space group != CIF-cell space group
    # group theory (Gamma point)
    n_vibrational_modes: int = 0
    n_vib_distinct: int = 0
    n_raman_active_modes: int = 0
    n_ir_active_modes: int = 0
    n_raman_ir_overlap: int = 0
    n_silent_modes: int = 0
    n_raman_dof: int = 0
    n_ir_dof: int = 0
    raman_mode_fraction: float = _NAN
    degeneracy_fraction: float = _NAN
    n_raman_irreps: int = 0
    max_raman_irrep_mult: int = 0
    raman_channel_entropy: float = _NAN
    group_theory_status: str = "not_run"
    # chemistry
    sigma_chi: float = _NAN
    mass_contrast: float = _NAN
    # bonds / local structure
    n_bonds: int = 0
    bonds_per_atom: float = _NAN
    mean_bond: float = _NAN
    cv_bond: float = _NAN               # std/mean of first-shell bond lengths, site-averaged
    bond_distortion: float = _NAN       # mean|d-<d>|/<d>, site-averaged (this was called cv_bond before)
    mean_cn: float = _NAN
    cn_std: float = _NAN
    mean_delta_chi_bond: float = _NAN
    mean_reduced_mass: float = _NAN
    n_anion_group_types: int = 0
    anion_group_o_fraction: float = _NAN
    polymerization_index: float = _NAN
    # SOAP
    soap_intra_species_dispersion: float = _NAN
    soap_n_env: int = 0
    soap_min_env_sep: float = _NAN
    # status
    ok: bool = False
    error: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


def _err(feat: StructureFeatures, msg: str) -> None:
    feat.error = f"{feat.error}; {msg}" if feat.error else msg


# --------------------------------------------------------------------------------------
# Main entry
# --------------------------------------------------------------------------------------
def compute_structure_features(cif_path: Path) -> StructureFeatures:
    cif_path = Path(cif_path)
    stem = cif_path.stem
    parts = stem.rsplit("__", 1)
    feat = StructureFeatures(mineral=parts[0] if len(parts) == 2 else stem,
                             amcsd_id=parts[1] if len(parts) == 2 else "",
                             cif_path=str(cif_path), formula="")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            _compute(cif_path, feat)
        except Exception as e:                                   # never let one CIF kill a batch
            _err(feat, f"unexpected: {type(e).__name__}: {e}")
    return feat


def _compute(cif_path: Path, feat: StructureFeatures) -> None:
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    try:
        structure, flags = load_structure(cif_path)
    except Exception as e:
        _err(feat, f"parse_failed: {e}")
        return
    feat.cif_repaired, feat.n_unknown_labels = flags["repaired"], flags["n_unknown_labels"]
    if feat.n_unknown_labels:
        _err(feat, f"{feat.n_unknown_labels} atom labels replaced by dummy X")
    n_fu = 0.0
    try:
        n_fu = formula_units(structure)
    except Exception as e:
        _err(feat, f"formula_units_failed: {e}")
    if structure.is_ordered or n_fu <= 0:
        feat.formula = structure.composition.reduced_formula
    else:                                   # per-formula-unit, occupancy-weighted, e.g. Mg1.8Fe0.2SiO4
        feat.formula = (structure.composition / n_fu).formula.replace(" ", "")
    feat.n_sites_conventional = len(structure)
    for k, v in disorder_features(structure).items():
        setattr(feat, k, v)

    # ---- symmetry -------------------------------------------------------------------
    prim = structure
    prim_is_idealised = False
    try:
        sga = SpacegroupAnalyzer(structure, symprec=SYMPREC, angle_tolerance=ANGLE_TOLERANCE)
        feat.space_group_symbol = sga.get_space_group_symbol()
        feat.space_group_number = int(sga.get_space_group_number())
        feat.crystal_system = sga.get_crystal_system()
        feat.point_group = sga.get_point_group_symbol()
        try:
            feat.n_wyckoff_orbits = len(sga.get_symmetrized_structure().equivalent_indices)
        except Exception as e:
            _err(feat, f"wyckoff_failed: {e}")
        try:
            prim = sga.get_primitive_standard_structure()
            prim_is_idealised = True
        except Exception:
            prim = structure.get_primitive_structure()
    except Exception as e:
        _err(feat, f"symmetry_failed: {type(e).__name__}: {e}")
        try:
            prim = structure.get_primitive_structure()
        except Exception:
            prim = structure
    feat.n_atom_primitive = len(prim)
    feat.z_primitive = n_fu * len(prim) / len(structure) if n_fu > 0 and len(structure) else _NAN

    # ---- group theory on the primitive cell --------------------------------------------
    if feat.space_group_number:
        try:
            prim_sga = SpacegroupAnalyzer(prim, symprec=SYMPREC, angle_tolerance=ANGLE_TOLERANCE)
            if prim_sga.get_space_group_number() != feat.space_group_number:
                feat.symmetry_consistent = False
                _err(feat, f"sg_mismatch_after_primitive({prim_sga.get_space_group_number()})")
            gt = compute_group_theory_features(
                prim, feat.point_group, prim_sga.get_symmetry_operations(cartesian=False))
        except Exception as e:
            gt = GroupTheoryFeatures(status="failed", error=f"{type(e).__name__}: {e}")
        feat.group_theory_status = gt.status
        if gt.error:
            _err(feat, f"group_theory_failed: {gt.error}")
        for k in ("n_vibrational_modes", "n_vib_distinct", "n_raman_active_modes", "n_ir_active_modes",
                  "n_raman_ir_overlap", "n_silent_modes", "n_raman_dof", "n_ir_dof",
                  "raman_mode_fraction", "degeneracy_fraction", "n_raman_irreps",
                  "max_raman_irrep_mult", "raman_channel_entropy", "is_centrosymmetric"):
            setattr(feat, k, getattr(gt, k))
    else:
        feat.group_theory_status = "skipped_no_symmetry"

    # ---- composition ---------------------------------------------------------------------
    try:
        comp = structure.composition.fractional_composition.get_el_amt_dict()
        syms = list(comp)
        fr = np.array([comp[s] for s in syms], float)
        chi = np.array([_element_data(s)[0] for s in syms])
        m = np.array([_element_data(s)[1] for s in syms])
        fr = fr / fr.sum()
        feat.sigma_chi = float(np.sqrt(np.sum(fr * (chi - np.sum(fr * chi)) ** 2)))
        feat.mass_contrast = float(np.log(m.max() / m.min()))
        feat.family = classify_family(structure.composition)
    except Exception as e:
        _err(feat, f"chem_failed: {e}")

    # ---- bonds ---------------------------------------------------------------------------
    try:
        for k, v in bond_features(structure).items():
            setattr(feat, k, v)
    except Exception as e:
        _err(feat, f"bond_failed: {type(e).__name__}: {e}")

    # ---- SOAP ----------------------------------------------------------------------------
    try:
        s = soap_features(prim)
        feat.soap_intra_species_dispersion = s["intra_species_dispersion"]
        feat.soap_n_env = s["n_env"]
        feat.soap_min_env_sep = s["min_env_sep"]
        if not prim_is_idealised:
            _err(feat, "soap_on_non_idealised_cell")
    except Exception as e:
        _err(feat, f"soap_failed({type(e).__name__}): {str(e)[:150]}")

    # 'ok' = the symmetry / group-theory core succeeded (other blocks report via `error`)
    feat.ok = bool(feat.space_group_number) and feat.group_theory_status == "ok"


# --------------------------------------------------------------------------------------
# Batch helpers / CLI
# --------------------------------------------------------------------------------------
def iter_cif_files(cif_dir: Path):
    yield from sorted(Path(cif_dir).glob("*.cif"))


def extract_directory(cif_dir: Path, out_csv: Path, n_jobs: int = 1) -> int:
    files = list(iter_cif_files(cif_dir))
    if n_jobs > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=n_jobs) as ex:
            feats = list(ex.map(compute_structure_features, files, chunksize=8))
    else:
        feats = [compute_structure_features(f) for f in files]
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=[f.name for f in fields(StructureFeatures)])
        w.writeheader()
        for f in feats:
            w.writerow(f.to_dict())
    return len(feats)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Extract structural features from CIF files")
    ap.add_argument("cif_dir", type=Path)
    ap.add_argument("-o", "--out", type=Path, default=Path("structure_features.csv"))
    ap.add_argument("-j", "--jobs", type=int, default=1)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO)
    log.info("wrote %d rows to %s", extract_directory(a.cif_dir, a.out, a.jobs), a.out)