"""
Raman spectral feature extraction from RRUFF processed Raman-spectrum text
files.

For each spectrum this computes the observable Raman-side descriptors used
as response variables:

    N_peak        -- number of resolved peaks
    S_peak        -- Shannon entropy of the peak-area distribution
    Gamma         -- median width at half prominence of resolved peaks
    dOmega_median -- median spacing between adjacent peak positions
    R_overlap     -- Gamma / dOmega_median  (band-overlap / resolvability index)
    shape_w1_distance -- 1-D Wasserstein (earth-mover's) distance between the
                         background-corrected spectral-intensity distribution
                         and a fixed reference distribution (see
                         `w1_distance_to_uniform_reference` for the exact
                         definition). Larger values indicate spectral weight
                         concentrated away from a flat/uniform profile.
    shape_centroid,
    shape_skewness    -- intensity-weighted centre of mass and asymmetry of the
                         SAME observed distribution used for shape_w1_distance.
                         These are peak-detection-free, so they remain
                         well-defined even for spectra with weak or ambiguous
                         discrete peaks.

Analysis window: peaks are only searched for in the 100-1300 cm^-1 "finger-
print" region.

Reference distribution for shape_w1_distance
---------------------------------------
`shape_w1_distance` is the 1-D Wasserstein distance between two distributions over
wavenumber, both restricted to the analysis window [WINDOW_MIN, WINDOW_MAX]
= [100, 1300] cm^-1:

  (i)  the OBSERVED distribution: the wavenumber grid points `x` retained
       after windowing, each weighted by the background-corrected intensity
       `y_bc` at that point (rolling-minimum baseline subtracted, clipped at
       zero; see `_read_spectrum` / baseline step above). This is the
       spectrum's own intensity profile treated as an (unnormalised)
       measure over wavenumber.
  (ii) the REFERENCE distribution: `n` points spaced evenly across the full
       nominal window [WINDOW_MIN, WINDOW_MAX], where `n = len(y_bc)` (i.e.
       the same number of support points as the observed spectrum), each
       given EQUAL weight (scipy's `wasserstein_distance` default when no
       `v_weights` is supplied). In the limit of dense sampling this
       converges to the continuous Uniform(WINDOW_MIN, WINDOW_MAX)
       distribution, i.e. a perfectly flat spectral-weight profile spanning
       the whole fingerprint region.

Two points are worth stating explicitly for reproducibility:
  * The reference always spans the FULL nominal window [100, 1300] cm^-1,
    regardless of the wavenumber range actually populated by a given
    spectrum's non-zero intensity -- it is not renormalised per spectrum.
  * `shape_w1_distance` is therefore best read as a measure of how unevenly
    spectral intensity is distributed across the fingerprint window
    relative to a flat baseline, not as a distance to any other measured
    or simulated spectrum.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks, peak_widths, savgol_filter
from scipy.stats import wasserstein_distance

WINDOW_MIN = 100.0
WINDOW_MAX = 1300.0
MIN_HEIGHT_FRAC = 0.02
MIN_PROMINENCE_FRAC = 0.015
MIN_PEAK_DISTANCE_CM = 3.0
N_SIGMA_NOISE = 5.0

# 定义物理频段
BANDS = {"low": (100, 400), "mid": (400, 700), "high": (700, 1300)}

@dataclass
class RamanFeatures:
    mineral: str
    rruff_id: str
    file_path: str
    quality: str
    orientation: str
    wavelength: str
    n_points: int = 0
    n_peak: int = 0
    s_peak: float = float("nan")
    gamma: float = float("nan")
    domega_median: float = float("nan")
    r_overlap: float = float("nan")
    n_low: int = 0
    s_low: float = float("nan")
    gamma_low: float = float("nan")
    n_mid: int = 0
    s_mid: float = float("nan")
    gamma_mid: float = float("nan")
    n_high: int = 0
    s_high: float = float("nan")
    gamma_high: float = float("nan")
    shape_w1_distance: float = float("nan")  # see w1_distance_to_uniform_reference() below for definition
    shape_centroid: float = float("nan")  # see spectral_moment_features() below for definition
    shape_skewness: float = float("nan")
    peak_positions_cm1: str = ""
    snr: float = float("nan")
    ok: bool = False
    error: str = ""

@dataclass(frozen=True)
class RamanSettings:
    baseline_window_fraction: float = 1 / 40
    smoothing_window_fraction: float = 1 / 100
    min_height_frac: float = MIN_HEIGHT_FRAC
    min_prominence_frac: float = MIN_PROMINENCE_FRAC
    min_peak_distance_cm: float = MIN_PEAK_DISTANCE_CM
    n_sigma_noise: float = N_SIGMA_NOISE


def w1_distance_to_uniform_reference(x: np.ndarray, y_bc: np.ndarray,
                                      window_min: float = WINDOW_MIN,
                                      window_max: float = WINDOW_MAX) -> float:
    """1-D Wasserstein distance from the observed intensity profile to a
    uniform reference spanning the analysis window.

    Parameters
    ----------
    x : wavenumber grid points (cm^-1) retained after windowing to
        [window_min, window_max]; must be strictly increasing.
    y_bc : background-corrected intensity at each point in `x` (same
        length as `x`); used as unnormalised weights, so absolute scale
        does not matter (scipy normalises internally).
    window_min, window_max : bounds of the analysis / reference window,
        defaults to the module-level fingerprint window [100, 1300] cm^-1.

    Reference distribution
    -----------------------
    The reference is `len(y_bc)` points spaced evenly across
    [window_min, window_max], each with equal weight -- i.e. a discretised
    Uniform(window_min, window_max) distribution using as many support
    points as the observed spectrum. It always spans the full nominal
    window regardless of where a given spectrum's intensity happens to be
    concentrated, so `shape_w1_distance` measures departure from a flat
    spectral-weight profile across the whole fingerprint region, not a
    distance to another measured spectrum.

    Returns
    -------
    float : the Wasserstein-1 distance, in wavenumber units (cm^-1).
    """
    reference_points = np.linspace(window_min, window_max, len(y_bc))
    return float(wasserstein_distance(x, reference_points, u_weights=y_bc))


def spectral_moment_features(x: np.ndarray, y_bc: np.ndarray) -> dict:
    """Weighted central-moment shape descriptors of the spectrum, computed
    on exactly the same (x, y_bc) observed distribution used by
    `w1_distance_to_uniform_reference` -- i.e. the background-corrected
    intensity at each retained wavenumber point, treated as an unnormalised
    measure over wavenumber within the fingerprint window. Together with
    `shape_w1_distance`, these give a peak-detection-free, model-free summary of
    spectral shape: where the spectral weight is centred, how broadly it is
    and whether it is asymmetric.

    Parameters
    ----------
    x : wavenumber grid points (cm^-1) retained after windowing.
    y_bc : background-corrected intensity at each point in `x` (same
        length as `x`); used as unnormalised weights.

    Returns
    -------
    dict with keys:
        shape_centroid : intensity-weighted mean wavenumber (cm^-1) -- the
            "centre of mass" of the spectrum within the fingerprint window.
        shape_skewness : standardised third central moment (Fisher-Pearson
            skewness). Positive values indicate a tail of spectral weight
            extending toward higher wavenumber than the bulk of the
            intensity, negative values toward lower wavenumber.

    `shape_skewness` is set to NaN (rather than 0 or an arbitrary fallback)
    when the weighted distribution is numerically degenerate, since the standardised
    moment is undefined in that limit.
    """
    x = np.asarray(x, dtype=float)
    w = np.asarray(y_bc, dtype=float)
    total = float(w.sum())
    nan = float("nan")
    if total <= 0 or not np.isfinite(total):
        return dict(shape_centroid=nan, shape_skewness=nan)

    p = w / total
    mu = float(np.sum(p * x))
    var = float(np.sum(p * (x - mu) ** 2))
    std = var ** 0.5
    if not np.isfinite(std) or std <= 1e-8:
        return dict(shape_centroid=mu, shape_skewness=nan)

    m3 = float(np.sum(p * (x - mu) ** 3))
    skewness = m3 / std ** 3
    return dict(shape_centroid=mu, shape_skewness=skewness)


def compute_band_features(positions: np.ndarray, fwhm: np.ndarray, areas: np.ndarray,
                           p_min: float, p_max: float) -> dict:
    """Per-band summary of detected peaks in [p_min, p_max), except the final band includes its upper bound.

    Parameters
    ----------
    positions, fwhm, areas : per-peak arrays (same length, same order) as computed
        in `compute_raman_features` -- peak wavenumber, half-prominence width, and area.
    p_min, p_max : band bounds in cm^-1 (see the module-level `BANDS` dict).

    Returns
    -------
    dict with keys "n" (peak count in the band), "s" (Shannon entropy of the
    band's peak-area distribution), "gamma" (median half-prominence width).
    Returns {"n": 0, "s": 0.0, "gamma": 0.0} when no peaks fall in the band.
    """
    mask = ((positions >= p_min) &
            ((positions <= p_max) if p_max == WINDOW_MAX else (positions < p_max)))
    if mask.sum() == 0:
        return {"n": 0, "s": 0.0, "gamma": 0.0}
    b_pos, b_fwhm, b_areas = positions[mask], fwhm[mask], areas[mask]
    p = b_areas / b_areas.sum()
    return {
        "n": len(b_pos),
        "s": float(-np.sum(p * np.log(p))),
        "gamma": float(np.median(b_fwhm))
    }


def _parse_meta_from_filename(path: Path):
    parts = path.stem.split("__")
    mineral = parts[0] if parts else path.stem
    rruff_id = parts[1] if len(parts) > 1 else ""
    wavelength = ""
    if len(parts) > 3:
        wavelength = parts[3]
    return mineral, rruff_id, wavelength

def _read_spectrum(path: Path):
    xs, ys = [], []
    with open(path, "r", encoding="latin-1") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("##"):
                continue
            m = re.match(r"^\s*([-\d.Ee+]+)\s*,\s*([-\d.Ee+]+)\s*$", line)
            if not m:
                continue
            try:
                x, y = float(m.group(1)), float(m.group(2))
            except ValueError:
                continue
            xs.append(x)
            ys.append(y)
    return np.array(xs), np.array(ys)

def compute_raman_features(path: Path, quality: str, orientation: str,
                           settings: RamanSettings | None = None) -> RamanFeatures:
    settings = settings or RamanSettings()
    mineral, rruff_id, wavelength = _parse_meta_from_filename(path)
    feat = RamanFeatures(
        mineral=mineral, rruff_id=rruff_id, file_path=str(path),
        quality=quality, orientation=orientation, wavelength=wavelength,
    )
    try:
        x, y = _read_spectrum(path)
        if len(x) < 20:
            feat.error = "too_few_points"
            return feat
        order = np.argsort(x)
        x, y = x[order], y[order]
        mask = (x >= WINDOW_MIN) & (x <= WINDOW_MAX)
        if mask.sum() < 20:
            feat.error = "window_too_short"
            return feat
        x, y = x[mask], y[mask]
        feat.n_points = len(x)

        dx_values = np.diff(x)
        dx_values = dx_values[dx_values > 0]
        if len(dx_values) == 0:
            feat.error = "invalid_wavenumber_spacing"
            return feat

        win = max(11, int(len(y) * settings.baseline_window_fraction) | 1)
        rolling_min = np.array([
            y[max(0, i - win): i + win + 1].min() for i in range(len(y))
        ])
        y_bc = np.clip(y - rolling_min, 0, None)

        dx = float(np.median(dx_values))
        if not np.isfinite(dx) or dx <= 0:
            feat.error = "invalid_wavenumber_spacing"
            return feat
        max_signal = y_bc.max()
        if max_signal <= 0:
            feat.error = "flat_spectrum"
            return feat

        sg_window = max(5, int(len(y_bc) * settings.smoothing_window_fraction) | 1)
        sg_window = min(sg_window, len(y_bc) - 1 if len(y_bc) % 2 == 0 else len(y_bc))
        if sg_window < 5:
            sg_window = 5
        if sg_window % 2 == 0:
            sg_window += 1
        try:
            y_smooth = savgol_filter(y_bc, window_length=sg_window, polyorder=3)
        except Exception:
            y_smooth = y_bc
        y_smooth = np.clip(y_smooth, 0, None)
        noise = y_bc - y_smooth
        sigma_noise = 1.4826 * float(np.median(np.abs(noise - np.median(noise))))
        feat.snr = float(max_signal / sigma_noise) if sigma_noise > 0 else float("inf")

        min_distance_pts = max(1, int(round(settings.min_peak_distance_cm / dx)))
        prominence_floor = max(settings.min_prominence_frac * max_signal, settings.n_sigma_noise * sigma_noise)
        height_floor = max(settings.min_height_frac * max_signal, settings.n_sigma_noise * sigma_noise)
        peaks, _ = find_peaks(
            y_smooth,
            height=height_floor,
            prominence=prominence_floor,
            distance=min_distance_pts,
        )
        if len(peaks) == 0:
            feat.error = "no_peaks"
            return feat

        widths_res = peak_widths(y_smooth, peaks, rel_height=0.5)
        widths_pts, _, left_ips, right_ips = widths_res
        fwhm = widths_pts * dx

        wide_res = peak_widths(y_smooth, peaks, rel_height=0.95)
        _, _, left95, right95 = wide_res
        
        trapz_func = getattr(np, "trapezoid", getattr(np, "trapz", None))
        areas = []
        for k in range(len(peaks)):
            lo = int(np.floor(left95[k]))
            hi = int(np.ceil(right95[k]))
            lo = max(0, lo)
            hi = min(len(x) - 1, hi)
            if hi <= lo:
                areas.append(float(y_bc[peaks[k]] * dx))
            else:
                if trapz_func is not None:
                    area_val = float(trapz_func(y_bc[lo:hi + 1], x[lo:hi + 1]))
                else:
                    area_val = float(np.sum(y_bc[lo:hi + 1]) * dx)
                areas.append(area_val)
        areas = np.array(areas)
        areas = np.clip(areas, 1e-12, None)

        positions = x[peaks]
        order2 = np.argsort(positions)
        positions = positions[order2]
        fwhm = fwhm[order2]
        areas = areas[order2]

        feat.n_peak = len(positions)
        feat.peak_positions_cm1 = ";".join(f"{position:.2f}" for position in positions)
        p = areas / areas.sum()
        feat.s_peak = float(-np.sum(p * np.log(p)))
        feat.gamma = float(np.median(fwhm))
        if len(positions) >= 2:
            feat.domega_median = float(np.median(np.diff(positions)))
            if feat.domega_median > 0:
                feat.r_overlap = feat.gamma / feat.domega_median
                
        for band_name, (p_min, p_max) in BANDS.items():
            b_feats = compute_band_features(positions, fwhm, areas, p_min, p_max)
            setattr(feat, f"n_{band_name}", b_feats["n"])
            setattr(feat, f"s_{band_name}", b_feats["s"])
            setattr(feat, f"gamma_{band_name}", b_feats["gamma"])

        feat.shape_w1_distance = w1_distance_to_uniform_reference(x, y_bc)
        for _k, _v in spectral_moment_features(x, y_bc).items():
            setattr(feat, _k, _v)

        feat.ok = True
    except Exception as e:
        feat.error = f"failed: {e}"
    return feat