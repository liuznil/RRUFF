"""
Raman spectral feature extraction from RRUFF processed Raman-spectrum text
files.

For each spectrum this computes the observable Raman-side descriptors used
as response variables:

    N_peak        -- number of resolved peaks
    S_peak        -- Shannon entropy of the peak-area distribution
    Gamma         -- median FWHM of resolved peaks (linewidth)
    dOmega_median -- median spacing between adjacent peak positions
    R_overlap     -- Gamma / dOmega_median  (band-overlap / resolvability index)
    w1_distance   -- Wasserstein distance to uniform distribution

Analysis window: peaks are only searched for in the 100-1300 cm^-1 "finger-
print" region.
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
    w1_distance: float = float("nan")
    snr: float = float("nan")
    ok: bool = False
    error: str = ""

def compute_band_features(positions, fwhm, areas, p_min, p_max):
    mask = (positions >= p_min) & (positions < p_max)
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

def compute_raman_features(path: Path, quality: str, orientation: str) -> RamanFeatures:
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

        win = max(11, len(y) // 40 | 1)
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

        sg_window = max(5, (len(y_bc) // 100) | 1)
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

        min_distance_pts = max(1, int(round(MIN_PEAK_DISTANCE_CM / dx)))
        prominence_floor = max(MIN_PROMINENCE_FRAC * max_signal, N_SIGMA_NOISE * sigma_noise)
        height_floor = max(MIN_HEIGHT_FRAC * max_signal, N_SIGMA_NOISE * sigma_noise)
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

        ideal_uniform = np.linspace(WINDOW_MIN, WINDOW_MAX, len(y_bc))
        feat.w1_distance = float(wasserstein_distance(x, ideal_uniform, u_weights=y_bc))

        feat.ok = True
    except Exception as e:
        feat.error = f"failed: {e}"
    return feat