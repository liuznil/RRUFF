#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Structural and composition-derived descriptors  ->  Raman-spectrum features : publication-grade statistical workflow
==========================================================================================
整合并升级  (严谨的回归 / 群论掩码 / LOFO) 与 (质控漏斗 / 偏相关 / 置换重要性 / Lasso / PLS) .

输入 : results/dataset.csv   (每行一种矿物: RRUFF AMCSD 结构/组成描述符 + RRUFF Raman 拉曼谱图特征)
输出 : <out>/figures/main/*.png (7 main figures), <out>/figures/supplementary/*.png (6 supplementary figures),
       <out>/tables/*.csv, <out>/results_summary.md,
       <out>/methods_text.md (方法段落, 数值由本次运行自动填入), <out>/run_config.json
依赖 : numpy pandas scipy scikit-learn matplotlib

分析链条 (每一步对应一个论文小节)
  S0  质控与队列描述        质控漏斗 + Table 1 + 化学族谱画像                                  -> Fig1, Table1
  S1  单变量关联            Spearman + 偏相关(控制 ln N_atom, ln SNR), Fisher 95%CI, BH-FDR  -> Fig2
  S2  分类因素              Kruskal-Wallis 效应量 eta^2_H (化学族/晶系/对称性/波长)          -> Fig3
  S3  代表性关系            每个谱特征自动选偏相关最强的结构变量, 散点 + 分箱中位线          -> Fig4
  S4  群论约束              N_Raman vs N_peak，检验群论模式数对实验峰数的可解释上限             -> Fig5
  S5  多元回归              M1 core / M2 extended / M3 extended+化学族固定效应;
                            z-score 系数, 按化学式聚类稳健 SE, VIF, BH-FDR                   -> Fig6
  S6  预测与增量价值        按化学式分组的重复 K 折 CV + cluster bootstrap；与测量/类别基线比较    -> Fig7, FigS1
  S7  重要性                留出折置换重要性 (单特征 / 特征组联合) + 特征组消融 dR2           -> Fig8
  S8  稀疏关系              嵌套分组 CV 的 Lasso(1-SE) + bootstrap 稳定性选择                -> Fig9
  S9  潜变量/稳健性         PLS + QC sensitivity；主要结果作为补充验证而非核心证据             -> Fig10, Table S
  S10 稳健性                在不同 SNR / 最少峰数阈值下重复偏相关, 报告符号/显著性一致率      -> Figure S6 + Table S


用法
  python data_analysis.py --csv results/dataset.csv --out results
  python data_analysis.py --fast          # 调试: 1 次重复 CV, 少量 bootstrap, 跳过消融/LOFO
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
import warnings
from collections import namedtuple
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import warnings
warnings.filterwarnings(
    "ignore",
    message=r".*sklearn\.utils\.parallel\.delayed.*"
)
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle, Ellipse, Polygon
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from scipy import stats  # noqa: E402
from scipy.stats import rankdata  # noqa: E402
from sklearn.cross_decomposition import PLSRegression  # noqa: E402
from sklearn.ensemble import HistGradientBoostingRegressor  # noqa: E402
from sklearn.linear_model import Lasso, lasso_path  # noqa: E402

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# =============================================================================
# 0. 配置: 特征 / 谱特征 / 分组 / 配色
# =============================================================================
GROUPS = ["size", "sym", "bond", "chem", "soap"]
GROUP_LABEL = {"size": "Size / complexity", "sym": "Symmetry / group theory", "bond": "Bond geometry",
               "chem": "Composition-derived chemistry", "soap": "Local environment (SOAP)"}
# Okabe-Ito 色盲友好配色
GROUP_COLOR = {"size": "#8C8C8C", "sym": "#0072B2", "bond": "#D55E00", "chem": "#009E73", "soap": "#CC79A7"}
FAM_PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#8C564B",
               "#7F7F7F", "#882255", "#332288", "#999933", "#000000"]

Feat = namedtuple("Feat", "name label group tf")     # tf: None | 'log' (ln, 非正->NaN) | 'log1' (ln(1+x))
FEATURES = [
    Feat("n_atom_primitive", "N_atom (primitive)", "size", "log"),
    Feat("n_wyckoff_orbits", "Wyckoff orbits", "size", "log1"),
    Feat("is_centrosymmetric", "Centrosymmetric (0/1)", "sym", None),
    Feat("n_raman_active_modes", "N_Raman (group theory)", "sym", "log"),   # 0 -> NaN, 不伪造 ln1=0
    Feat("raman_mode_fraction", "Raman-active fraction", "sym", None),
    Feat("degeneracy_fraction", "Degeneracy fraction", "sym", None),
    Feat("raman_channel_entropy", "Raman irrep entropy", "sym", None),
    Feat("n_raman_irreps", "N Raman irreps", "sym", None),
    Feat("max_raman_irrep_mult", "Max Raman irrep mult.", "sym", None),
    Feat("mean_bond", "Mean bond length", "bond", None),
    Feat("cv_bond", "Bond-length CV", "bond", None),
    Feat("bond_distortion", "Bond distortion", "bond", None),
    Feat("mean_cn", "Mean coordination no.", "bond", None),
    Feat("cn_std", "CN std", "bond", None),
    Feat("bonds_per_atom", "Bonds per atom", "bond", None),
    Feat("mean_delta_chi_bond", "Bond electroneg. diff.", "bond", None),
    Feat("mean_reduced_mass", "Bond reduced mass", "bond", "log"),
    Feat("sigma_chi", "Electroneg. dispersion", "chem", None),
    Feat("mass_contrast", "Mass contrast", "chem", None),
    Feat("anion_group_o_fraction", "Oxyanion fraction", "chem", None),
    Feat("n_anion_group_types", "N anion-group types", "chem", None),
    Feat("polymerization_index", "Polyhedral polymerisation", "chem", None),
    Feat("mixing_entropy", "Site-mixing entropy", "chem", None),
    Feat("soap_intra_species_dispersion", "SOAP intra-species disp.", "soap", None),
    Feat("soap_n_env", "SOAP N environments", "soap", "log1"),
    Feat("soap_min_env_sep", "SOAP min env. separation", "soap", "log"),
]
FEAT = {f.name: f for f in FEATURES}

Targ = namedtuple("Targ", "name label tf")            # tf 额外允许 'auto' (按偏度数据驱动决定是否取 ln)
TARGETS = [
    Targ("n_peak", "N_peak", "log"),
    Targ("evenness", "Peak-intensity evenness", None),
    Targ("gamma", "Mean linewidth gamma", "log"),
    Targ("domega_median", "Median peak spacing", "log"),
    Targ("r_overlap", "Peak overlap ratio", "auto"),
    Targ("frac_low", "Low-band fraction", None),
    Targ("frac_mid", "Mid-band fraction", None),
    Targ("frac_high", "High-band fraction", None),
    Targ("w1_distance", "W1 distance", "log"),
    Targ("eta", "Mapping efficiency eta", None),
]
TARG = {t.name: t for t in TARGETS}
PLS_EXCLUDE = {"eta"}
# eta = N_peak / N_Raman 是 N_Raman 的比值: 二者相关是定义性的, 不做检验; eta 也不进入 ML 环节 (特征会直接泄漏答案)
TRIVIAL_PAIRS = {("n_raman_active_modes", "eta")}
ML_EXCLUDE = {"eta"}

# 回归预测变量集合 (需近似完整, 不含结构性 NaN 的描述符)
CORE_OLS = ["n_atom_primitive", "cv_bond", "sigma_chi", "soap_intra_species_dispersion"]
EXT_OLS = CORE_OLS + ["degeneracy_fraction", "raman_channel_entropy", "mass_contrast",
                      "mean_cn", "mean_delta_chi_bond", "mean_reduced_mass"]
GT_RESPONSES = ["n_peak", "evenness"]                   # 群论子样本回归: N_Raman 代替 N_atom (近共线)

# 群论列: symmetry_consistent==False 时整体置 NaN (所描述的群不是所报告的群)
GROUP_THEORY_COLUMNS = ["n_vib_distinct", "n_raman_active_modes", "n_ir_active_modes", "n_raman_ir_overlap",
                        "n_silent_modes", "n_raman_dof", "n_ir_dof", "raman_mode_fraction", "degeneracy_fraction",
                        "n_raman_irreps", "max_raman_irrep_mult", "raman_channel_entropy"]
REQUIRED = ["n_atom_primitive", "family", "n_peak", "s_peak", "gamma", "r_overlap", "w1_distance"] + CORE_OLS[1:]

FAM_ABBR = {"phosphate_arsenate_vanadate": "phosphate/arsenate/vanadate", "sulfide_selenide_telluride": "sulfide/selenide/telluride",
            "oxide_hydroxide": "oxide/hydroxide", "native_element": "native element"}
FDR_ALPHA = 0.05


def fam_name(f: str) -> str:
    return FAM_ABBR.get(f, str(f).replace("_", " "))


def flabel(c: str) -> str:
    f = FEAT.get(c)
    return c if f is None else (("ln " if f.tf else "") + f.label)


def tlabel(c: str) -> str:
    t = TARG.get(c)
    return c if t is None else (("ln " if t.tf in ("log", "auto_log") else "") + t.label)


def fgroup(c):
    # 如果 c 已经是 Feat 对象，直接返回其 group 属性
    if isinstance(c, Feat):
        return c.group
    # 兼容旧的字符串键查找方式
    return FEAT[c].group

# =============================================================================
# 1. 数据读取 / 质控 / 特征构造
# =============================================================================
def _as_bool(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s
    m = s.astype(str).str.strip().str.lower().map({"true": True, "false": False, "1": True, "0": False,
                                                    "1.0": True, "0.0": False})
    return m.fillna(True).astype(bool)


def read_raw(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"{path} not found - run the dataset builder first")
    df = pd.read_csv(path).replace([np.inf, -np.inf], np.nan)
    if "soap_distortion" in df.columns and "soap_intra_species_dispersion" not in df.columns:
        raise SystemExit("dataset.csv was written by the OLD build (has 'soap_distortion'); rebuild it first.")
    miss = [c for c in REQUIRED if c not in df.columns]
    if miss:                                          # 核心列缺失直接报错, 不静默丢弃
        raise SystemExit(f"dataset.csv is missing required columns: {miss}")
    df["symmetry_consistent"] = _as_bool(df["symmetry_consistent"]) if "symmetry_consistent" in df else True
    gt = [c for c in GROUP_THEORY_COLUMNS if c in df.columns]
    df[gt] = df[gt].astype(float)
    df.loc[~df["symmetry_consistent"], gt] = np.nan
    if "is_centrosymmetric" in df.columns:
        df["is_centrosymmetric"] = _as_bool(df["is_centrosymmetric"]).astype(float)
    return df


def wl_group(w) -> str:
    if pd.isna(w):
        return "unknown"
    m = re.search(r"(\d+)", str(w))               # 兼容 "532 nm" / "532" / 532.0
    if not m:
        return "unknown"
    v = int(m.group(1))
    return "532 nm" if v == 532 else "780-785 nm" if 780 <= v <= 785 else "other"


def apply_qc(df: pd.DataFrame, qualities, min_snr, min_peaks, orientation, verbose=True):
    """质控漏斗. 缺失的质控列会被跳过并提示 (而不是报错)."""
    steps, d = [("All samples", len(df))], df

    def step(label, mask):
        nonlocal d
        d = d[mask(d)]
        steps.append((label, len(d)))

    if "raman_quality" in d and qualities:
        step(f"Spectrum quality in {{{','.join(qualities)}}}", lambda x: x["raman_quality"].isin(qualities))
    elif verbose:
        print("  [QC] raman_quality missing or not requested: skipped")
    if "raman_snr" in d and min_snr is not None:
        step(f"SNR >= {min_snr:g}", lambda x: x["raman_snr"] >= min_snr)
    elif verbose:
        print("  [QC] raman_snr missing: skipped")
    step(f"Detected peaks >= {min_peaks}", lambda x: x["n_peak"] >= min_peaks)
    if "raman_orientation" in d and orientation != "any":
        step(f"Orientation = {orientation}", lambda x: x["raman_orientation"] == orientation)
    return d.copy(), pd.DataFrame(steps, columns=["step", "n"])


def _auto_log(x: pd.Series) -> bool:
    x = x.dropna()
    if len(x) < 30 or x.min() <= 0:
        return False
    return abs(stats.skew(x)) > 1.0 and abs(stats.skew(np.log(x))) < abs(stats.skew(x))


def _tf(s: pd.Series, tf):
    if tf == "log":
        return np.log(s.where(s > 0))
    if tf == "log1":
        return np.log1p(s.clip(lower=0))
    return s


def build_frames(d: pd.DataFrame):
    """返回 (R 原始尺度含派生量, D 变换后, 变换记录)."""
    R = d.reset_index(drop=True).copy()
    with np.errstate(all="ignore"):
        R["evenness"] = R["s_peak"] / np.log(R["n_peak"].where(R["n_peak"] > 1))   # 去掉 s ≈ ln n 的平凡依赖
        for b in ("low", "mid", "high"):
            if f"n_{b}" in R:
                R[f"frac_{b}"] = R[f"n_{b}"] / R["n_peak"]
        if "n_raman_active_modes" in R:
            R["eta"] = R["n_peak"] / R["n_raman_active_modes"].where(R["n_raman_active_modes"] > 0)
    R = R.replace([np.inf, -np.inf], np.nan)
    R["wl_group"] = R["raman_wavelength"].map(wl_group) if "raman_wavelength" in R else "unknown"
    if "n_anion_group_types" in R:
        R["n_anion_group_types_cat"] = R["n_anion_group_types"].astype("Int64").astype(str)
    D, rec = R.copy(), []
    for f in FEATURES:
        if f.name in D:
            D[f.name] = _tf(D[f.name], f.tf)
            rec.append(("feature", f.name, f.tf or "none"))
    for t in TARGETS:
        if t.name in D:
            tf = t.tf
            if tf == "auto":
                tf = "log" if _auto_log(R[t.name]) else None
                TARG[t.name] = Targ(t.name, t.label, "auto_log" if tf else None)   # 让标签反映实际变换
            D[t.name] = _tf(D[t.name], tf)
            rec.append(("target", t.name, tf or "none"))
    if "raman_snr" in D:
        D["raman_snr"] = _tf(D["raman_snr"], "log")
    return R, D, pd.DataFrame(rec, columns=["kind", "variable", "transform"])


def group_column(D: pd.DataFrame) -> str:
    return "formula" if "formula" in D else ("mineral" if "mineral" in D else None)


# =============================================================================
# 2. 统计工具 (全部 numpy/scipy 实现)
# =============================================================================
def bh_fdr(p) -> np.ndarray:
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    pv, n = p[ok], ok.sum()
    if n == 0:
        return out
    o = np.argsort(pv)
    r = np.minimum.accumulate((pv[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    q = np.empty(n)
    q[o] = np.minimum(r, 1)
    out[ok] = q
    return out


def stars(q) -> str:
    return "" if q is None or not np.isfinite(q) else "***" if q < 1e-3 else "**" if q < 1e-2 else "*" if q < 5e-2 else ""


def fmt_p(p) -> str:
    return "" if not np.isfinite(p) else "<0.001" if p < 1e-3 else f"{p:.3f}"


def corr_table(D, feats, targets, covars=()):
    """Spearman 相关 (给 covars 则为秩偏相关). 成对删除缺失; Fisher-z (Bonett-Wright SE) 95%CI; BH-FDR."""
    covars, rows = list(covars), []
    for x in feats:
        for y in targets:
            nan_row = (x, y, np.nan, np.nan, np.nan, np.nan, 0)
            if x in covars or x == y or (x, y) in TRIVIAL_PAIRS:
                rows.append(nan_row)
                continue
            s = D[[x, y] + covars].dropna()
            n, k = len(s), len(covars)
            if n < 10 + k or s[x].nunique() < 2 or s[y].nunique() < 2:      # 二元变量(0/1)也允许
                rows.append((x, y, np.nan, np.nan, np.nan, np.nan, n))
                continue
            rx, ry = rankdata(s[x]), rankdata(s[y])
            if k:
                Z = np.column_stack([np.ones(n)] + [rankdata(s[c]) for c in covars])
                rx = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]
                ry = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
            r = np.corrcoef(rx, ry)[0, 1]
            dfree = n - 2 - k
            t = r * np.sqrt(dfree / max(1e-12, 1 - r ** 2))
            p = 2 * stats.t.sf(abs(t), dfree) if np.isfinite(r) else np.nan
            se = np.sqrt((1 + r ** 2 / 2) / max(1, n - 3 - k))
            z = np.arctanh(np.clip(r, -0.999999, 0.999999))
            rows.append((x, y, r, p, np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se), n))
    T = pd.DataFrame(rows, columns=["feature", "target", "rho", "p", "ci_lo", "ci_hi", "n"])
    T["q"] = bh_fdr(T["p"].values)
    return T


def to_mat(T, rows, cols, col):
    return T.pivot(index="feature", columns="target", values=col).reindex(index=rows, columns=cols).to_numpy(float)


def kruskal_eta2(D, factor, targets, min_n=8):
    """Kruskal-Wallis eta^2_H = (H - k + 1)/(n - k) (Tomczak & Tomczak 2014); 样本 < min_n 的水平剔除."""
    rows = []
    for y in targets:
        s = D[[factor, y]].dropna()
        vc = s[factor].value_counts()
        s = s[s[factor].isin(vc[vc >= min_n].index)]
        gs = [g[y].values for _, g in s.groupby(factor)]
        k, n = len(gs), len(s)
        if k < 2:
            rows.append((factor, y, np.nan, np.nan, n, k))
            continue
        H, p = stats.kruskal(*gs)
        rows.append((factor, y, max(0.0, (H - k + 1) / (n - k)), p, n, k))
    return pd.DataFrame(rows, columns=["factor", "target", "eta2_H", "p", "n", "k_levels"])


def ols_robust(X: np.ndarray, y: np.ndarray, names, groups=None):
    """OLS + 稳健 SE. groups 给定 -> 聚类稳健 CR1 (df=G-1), 否则 HC3 (df=n-k). X 已含常数列(第0列)."""
    n, k = X.shape
    XtXi = np.linalg.pinv(X.T @ X)
    beta = XtXi @ X.T @ y
    e = y - X @ beta
    G = None
    if groups is not None:
        codes = pd.factorize(groups)[0]
        G = codes.max() + 1
    if G is not None and G >= 20:
        S = np.zeros((G, k))
        np.add.at(S, codes, X * e[:, None])
        meat = S.T @ S
        cov = XtXi @ meat @ XtXi * (G / (G - 1)) * ((n - 1) / (n - k))
        dfree, kind = G - 1, "cluster"
    else:
        h = np.einsum("ij,jk,ik->i", X, XtXi, X)
        u = e / np.clip(1 - h, 1e-8, None)
        cov = XtXi @ (X.T * u ** 2) @ X @ XtXi
        dfree, kind = n - k, "HC3"
    se = np.sqrt(np.clip(np.diag(cov), 0, None))
    t = beta / se
    p = 2 * stats.t.sf(np.abs(t), dfree)
    tc = stats.t.ppf(0.975, dfree)
    sst = ((y - y.mean()) ** 2).sum()
    r2 = 1 - (e ** 2).sum() / sst if sst > 0 else np.nan
    return dict(names=list(names), beta=beta, se=se, p=p, lo=beta - tc * se, hi=beta + tc * se, r2=r2,
                adj_r2=1 - (1 - r2) * (n - 1) / (n - k), n=n, k=k, se_kind=kind, n_clusters=G, dfree=dfree)


def vif(Xz: np.ndarray) -> np.ndarray:
    C = np.corrcoef(Xz, rowvar=False)
    return np.diag(np.linalg.pinv(np.atleast_2d(C)))


def group_folds(groups, k, seed):
    """按 group 划分折 (同一化学式不跨折), 贪心均衡各折样本数. 返回 [(train_idx, test_idx)]."""
    groups = np.asarray(groups)
    sizes = pd.Series(groups).value_counts()
    k = int(min(k, len(sizes)))
    rng = np.random.default_rng(seed)
    order = sizes.index.to_numpy()[np.lexsort((rng.random(len(sizes)), -sizes.values))]
    load, fold_of = np.zeros(k), {}
    for g in order:
        i = int(load.argmin())
        fold_of[g] = i
        load[i] += sizes[g]
    f = np.array([fold_of[g] for g in groups])
    return [(np.where(f != i)[0], np.where(f == i)[0]) for i in range(k)]


def make_boot_idx(groups, n_boot, seed):
    """按 group (化学式) 聚类 bootstrap 的行索引, 可在多个模型间复用 (配对比较)."""
    idx_by_g = [np.asarray(v) for v in pd.Series(np.arange(len(groups))).groupby(np.asarray(groups)).apply(list)]
    rng = np.random.default_rng(seed)
    G = len(idx_by_g)
    return [np.concatenate([idx_by_g[i] for i in rng.integers(0, G, G)]) for _ in range(n_boot)]


def boot_r2(y, preds, boots):
    """preds: (R, n) 各重复的 OOF 预测. 返回 bootstrap 分布 (跨重复取平均 R2)."""
    preds = np.atleast_2d(preds)
    out = []
    for idx in boots:
        yy = y[idx]
        sst = ((yy - yy.mean()) ** 2).sum()
        if sst > 0:
            out.append(np.mean([1 - ((yy - p[idx]) ** 2).sum() / sst for p in preds]))
        else:
            out.append(np.nan)
    return np.asarray(out)


def r2_pooled(y, p):
    sst = ((y - y.mean()) ** 2).sum()
    return 1 - ((y - p) ** 2).sum() / sst if sst > 0 else np.nan


class Prep:
    """训练折内拟合的预处理: 中位数填补 + (缺失率>阈值的特征)缺失指示列(去重) + 标准化. 杜绝泄漏."""

    def __init__(self, na_thresh=0.02):
        self.na_thresh = na_thresh

    def fit(self, X: pd.DataFrame):
        self.cols = list(X.columns)
        self.med = X.median().fillna(0.0)
        miss = X.isna().mean()
        cand = [c for c in self.cols if miss[c] > self.na_thresh]
        self.ind = []
        if cand:
            ind = X[cand].isna().astype(float).T.drop_duplicates().T     # 同一缺失模式(如群论掩码)只保留一列
            self.ind = list(ind.columns)
        Z = self._raw(X)
        self.mu, self.sd = Z.mean(), Z.std(ddof=0).replace(0, 1.0)
        self.names = list(Z.columns)
        return self

    def _raw(self, X):
        Z = X[self.cols].fillna(self.med)
        if self.ind:
            Z = pd.concat([Z, X[self.ind].isna().astype(float).add_prefix("NA:")], axis=1)
        return Z

    def transform(self, X):
        return ((self._raw(X) - self.mu) / self.sd).to_numpy(float)


def make_hgb(seed, fast=False):
    return HistGradientBoostingRegressor(max_iter=80 if fast else 150, learning_rate=0.06, max_depth=3,
                                         min_samples_leaf=15, l2_regularization=1.0, early_stopping=False,
                                         random_state=seed)


# =============================================================================
# 3. 出图工具
# =============================================================================
def set_style():
    """Global publication style: 8-pt text, restrained lines."""
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8.0, "axes.labelsize": 8.0, "axes.titlesize": 8.5, "axes.titleweight": "bold",
        "xtick.labelsize": 7.2, "ytick.labelsize": 7.2, "legend.fontsize": 7.2,
        "axes.linewidth": 0.65, "lines.linewidth": 1.0, "patch.linewidth": 0.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 120, "savefig.dpi": 300, "savefig.bbox": "tight", "savefig.pad_inches": 0.03,
        "mathtext.default": "regular", "axes.unicode_minus": False})


class Saver:
    """Save publication figures as high-resolution PNG only."""

    def __init__(self, fig_dir: Path):
        self.dir = fig_dir
        fig_dir.mkdir(parents=True, exist_ok=True)

    def __call__(self, fig, name):
        fig.savefig(self.dir / f"{name}.png", dpi=DPI)
        plt.close(fig)


# 模块级开关 (main() 里由命令行覆盖)
FIG_TITLES = False          # 期刊图内通常不放总标题, 由图注承担; --fig-titles 可打开
DPI = 300                   # PNG 分辨率
QA_LOG: list = []           # 每张图的排版检查结果


def panel(ax, letter, dx=-0.06, dy=1.04):
    """登记面板字母. 有左对齐标题 -> 并入标题 "(a) Title" (不会重叠); 否则放在轴外左上角.
    字母在保存前(_finalize_figure)统一落实, 因此 panel() 与 set_title() 的调用先后无关."""
    ax._panel_letter, ax._panel_xy = letter, (dx, dy)
    if ax.get_title(loc="left"):
        _apply_panel(ax)


def _apply_panel(ax):
    letter = getattr(ax, "_panel_letter", None)
    if not letter or getattr(ax, "_panel_done", False):
        return
    t = ax.get_title(loc="left")
    if t:
        ax.set_title(f"({letter}) {t}", loc="left")
    else:
        dx, dy = ax._panel_xy
        ax.text(dx, dy, f"({letter})", transform=ax.transAxes, fontweight="bold", fontsize=9.5,
                va="bottom", ha="right", clip_on=False)
    ax._panel_done = True


def sig_marker(q):
    if q is None or not np.isfinite(q):
        return ""
    return "***" if q < 1e-3 else "**" if q < 1e-2 else "*" if q < 5e-2 else ""


def annotate_stats(ax, text, x=0.03, y=0.97, fontsize=7.2, ha="left"):
    ax.text(x, y, text, transform=ax.transAxes, va="top", ha=ha, fontsize=fontsize,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="0.78", alpha=0.92))


def fam_colors(D):
    fams = D["family"].value_counts().index.tolist()
    return {f: FAM_PALETTE[i % len(FAM_PALETTE)] for i, f in enumerate(fams)}


def group_legend(fig, loc="lower center", ncol=5, **kw):
    hs = [Patch(color=GROUP_COLOR[g], label=GROUP_LABEL[g]) for g in GROUPS]
    return fig.legend(handles=hs, loc=loc, ncol=ncol, frameon=False, **kw)


def feat_rows(feats):
    order = sorted(feats, key=lambda c: GROUPS.index(fgroup(c)))
    cols = [GROUP_COLOR[fgroup(c)] for c in order]
    seps = [i for i in range(1, len(order)) if fgroup(order[i]) != fgroup(order[i - 1])]
    return order, cols, seps


def draw_heat(ax, M, Q, xlabels, ylabels, ycolors=None, vmin=-0.6, vmax=0.6, cmap="RdBu_r", fmt="{:.2f}",
              title=None, show_y=True, annot_min=0.0, sep_rows=(), fs=6.3):
    cm = plt.get_cmap(cmap).copy()
    cm = cm.with_extremes(bad="#efefef")
    im = ax.imshow(np.ma.masked_invalid(M), cmap=cm, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(range(M.shape[1]), xlabels, rotation=40, ha="right")
    ax.set_yticks(range(M.shape[0]), ylabels if show_y else [""] * len(ylabels))
    if ycolors is not None and show_y:
        for t, c in zip(ax.get_yticklabels(), ycolors):
            t.set_color(c)
    lim = 0.62 * max(abs(vmin), abs(vmax))
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if not np.isfinite(v) or abs(v) < annot_min:
                continue
            s = fmt.format(v) + (stars(Q[i, j]) if Q is not None else "")
            ax.text(j, i, s, ha="center", va="center", fontsize=fs, color="white" if abs(v) > lim else "#222")
    for r in sep_rows:
        ax.axhline(r - 0.5, color="white", lw=2.2)
    ax.set_xticks(np.arange(-.5, M.shape[1]), minor=True)
    ax.set_yticks(np.arange(-.5, M.shape[0]), minor=True)
    ax.grid(which="minor", color="white", lw=0.6)
    ax.tick_params(which="both", length=0)
    for s_ in ax.spines.values():
        s_.set_visible(False)
    if title:
        ax.set_title(title, loc="left")
    return im


# =============================================================================
# 4. S0 队列描述
# =============================================================================
def iqr_str(s: pd.Series) -> str:
    s = s.dropna()
    return "" if s.empty else f"{s.median():.3g} [{s.quantile(.25):.3g}-{s.quantile(.75):.3g}]"


def step0_cohort(R, D, funnel, tg, save, tab):
    cols = [c for c in ["n_atom_primitive", "n_peak", "gamma", "domega_median", "r_overlap", "w1_distance", "eta"] if c in R]
    rows = [dict(group="All", n=len(R), **{c: iqr_str(R[c]) for c in cols})]
    for f, g in R.groupby("family"):
        rows.append(dict(group=fam_name(f), n=len(g), **{c: iqr_str(g[c]) for c in cols}))
    pd.DataFrame(rows).to_csv(tab / "table1_cohort_median_iqr.csv", index=False)
    funnel.to_csv(tab / "table_qc_funnel.csv", index=False)

    Z = (D[tg] - D[tg].mean()) / D[tg].std()
    Z["family"] = D["family"]
    prof = Z.groupby("family").median()
    cnt = D["family"].value_counts()
    prof = prof.loc[cnt.index]
    fig = plt.figure(figsize=(11.5, 4.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 2.0], wspace=0.9)
    ax = fig.add_subplot(gs[0])
    y = np.arange(len(funnel))[::-1]
    ax.barh(y, funnel["n"], color="#5b8db8")
    for yi, n in zip(y, funnel["n"]):
        ax.text(n + funnel["n"].max() * 0.015, yi, f"{n}", va="center", fontsize=7.5)
    ax.set_yticks(y, funnel["step"])
    ax.set_xlim(0, funnel["n"].max() * 1.15)
    ax.set_xlabel("Number of spectra")
    ax.set_title("Quality-control funnel", loc="left")
    panel(ax, "a", dx=-0.02)
    ax2 = fig.add_subplot(gs[1])
    im = draw_heat(ax2, prof.values, None, [tlabel(t) for t in tg], [f"{fam_name(f)} (n={cnt[f]})" for f in prof.index],
                   vmin=-1.2, vmax=1.2, title="Spectral profile by chemical family (median z-score)", fs=6.5)
    fig.colorbar(im, ax=ax2, fraction=0.03, pad=0.02).set_label("z-score")
    panel(ax2, "b", dx=-0.02)
    save(fig, "Fig1_cohort_and_family_profile")


# =============================================================================
# 5. S1 相关 / S2 分类 / S3 代表关系 / S4 群论
# =============================================================================
COVARS = ["n_atom_primitive", "raman_snr"]


def step1_correlation(D, feats, tg, save, tab, covars):
    raw = corr_table(D, feats, tg)
    par = corr_table(D, feats, tg, covars)
    raw.to_csv(tab / "table_corr_spearman.csv", index=False)
    par.to_csv(tab / "table_corr_partial.csv", index=False)
    order, cols, seps = feat_rows(feats)
    yl, xl = [flabel(c) for c in order], [tlabel(t) for t in tg]
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 8.6), gridspec_kw={"wspace": 0.05})
    for k, (ax, T, ttl) in enumerate(zip(axs, (raw, par), ("Spearman $\\rho$", "Partial $\\rho$ (adj. ln N_atom, ln SNR)"))):
        im = draw_heat(ax, to_mat(T, order, tg, "rho"), to_mat(T, order, tg, "q"), xl, yl, cols, title=ttl,
                       show_y=(k == 0), sep_rows=seps)
        panel(ax, "ab"[k], dx=-0.02 if k else -0.42)
    fig.colorbar(im, ax=axs, fraction=0.02, pad=0.02).set_label("correlation coefficient")
    fig.text(0.5, -0.035, "* q<0.05, ** q<0.01, *** q<0.001 (Benjamini-Hochberg over all cells of the panel); grey = covariate itself / undefined",
             ha="center", fontsize=7.5)
    group_legend(fig, loc="upper center", bbox_to_anchor=(0.5, 0.965))
    save(fig, "Fig2_correlation_raw_vs_partial")
    return raw, par


def step2_categorical(D, tg, save, tab):
    facs = [("family", "Chemical family"), ("crystal_system", "Crystal system"), ("is_centrosymmetric", "Centrosymmetry"),
            ("n_anion_group_types_cat", "N anion-group types"), ("wl_group", "Excitation wavelength (measurement)")]
    facs = [(f, n) for f, n in facs if f in D]
    res = pd.concat([kruskal_eta2(D, f, tg) for f, _ in facs], ignore_index=True)
    res["q"] = bh_fdr(res["p"].values)
    res.to_csv(tab / "table_categorical_effects_kruskal.csv", index=False)
    M = np.array([[res[(res.factor == f) & (res.target == t)]["eta2_H"].values[0] for t in tg] for f, _ in facs])
    Q = np.array([[res[(res.factor == f) & (res.target == t)]["q"].values[0] for t in tg] for f, _ in facs])
    pref = [t for t in ["n_peak", "gamma", "frac_high", "w1_distance"] if t in tg]
    fig = plt.figure(figsize=(11.5, 9.5))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.8, 1.1, 1.1], hspace=0.75, wspace=0.18)
    ax = fig.add_subplot(gs[0, :])
    im = draw_heat(ax, M, Q, [tlabel(t) for t in tg], [n for _, n in facs], vmin=0, vmax=max(0.3, np.nanmax(M)),
                   cmap="Purples", title="Variance in spectral features explained by categorical factors ($\\eta^2_H$)", fs=7)
    fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01).set_label("$\\eta^2_H$")
    panel(ax, "a", dx=-0.01)
    cf = fam_colors(D)
    for k, y in enumerate(pref):
        ax = fig.add_subplot(gs[1 + k // 2, k % 2])
        order_f = D.groupby("family")[y].median().sort_values().index.tolist()
        data = [D.loc[D.family == f, y].dropna().values for f in order_f]
        bp = ax.boxplot(data, patch_artist=True, showfliers=False, widths=0.6, medianprops=dict(color="k", lw=1.2))
        for patch, f in zip(bp["boxes"], order_f):
            patch.set_facecolor(cf[f])
            patch.set_alpha(0.85)
        ax.set_xticks(range(1, len(order_f) + 1), [fam_name(f) for f in order_f], rotation=40, ha="right")
        ax.set_ylabel(tlabel(y))
        ax.set_title(f"{tlabel(y)} by family", loc="left")
        panel(ax, "bcde"[k], dx=-0.09)
    save(fig, "Fig3_categorical_effects")
    return res


def step3_key_scatter(D, par, raw, feats, save, tab):
    show = [t for t in ["n_peak", "gamma", "domega_median", "frac_low", "frac_high", "evenness"] if t in D]
    used, picks = set(), []
    for y in show:
        sub = par[(par.target == y) & par.feature.isin(feats)].copy()
        sub = sub[sub.feature.map(lambda c: D[c].nunique() > 10)]
        sub["a"] = sub["rho"].abs()
        for _, r in sub.sort_values("a", ascending=False).iterrows():
            if r.feature not in used and np.isfinite(r.a):
                used.add(r.feature)
                picks.append((r.feature, y))
                break
    cf = fam_colors(D)
    nrow = int(np.ceil(len(picks) / 3))
    fig, axs = plt.subplots(nrow, 3, figsize=(11, 3.6 * nrow))
    rows = []
    for ax, (x, y), lab in zip(np.atleast_1d(axs).ravel(), picks, "abcdef"):
        s = D[[x, y, "family"]].dropna()
        ax.scatter(s[x], s[y], s=7, c=[cf[f] for f in s["family"]], alpha=0.5, lw=0)
        q = pd.qcut(s[x], 10, duplicates="drop")
        g = s.groupby(q, observed=True)
        xm = g[x].median()
        ax.fill_between(xm, g[y].quantile(.25), g[y].quantile(.75), color="k", alpha=0.13, lw=0)
        ax.plot(xm, g[y].median(), "k-", lw=1.8)
        rr = raw[(raw.feature == x) & (raw.target == y)].iloc[0]
        pp = par[(par.feature == x) & (par.target == y)].iloc[0]
        ax.text(0.03, 0.97, f"Spearman $\\rho$ = {rr.rho:+.2f}\npartial $\\rho$ = {pp.rho:+.2f} [{pp.ci_lo:+.2f}, {pp.ci_hi:+.2f}]\nn = {int(rr.n)}",
                transform=ax.transAxes, va="top", fontsize=7, bbox=dict(boxstyle="round", fc="white", ec="#bbb", alpha=0.9))
        ax.set_xlabel(flabel(x))
        ax.set_ylabel(tlabel(y))
        panel(ax, lab, dx=-0.13)
        rows.append((y, x, rr.rho, pp.rho, pp.ci_lo, pp.ci_hi, pp.q))
    for ax in np.atleast_1d(axs).ravel()[len(picks):]:
        ax.axis("off")
    cnt = D["family"].value_counts().index
    fig.legend(handles=[Patch(color=cf[f], label=fam_name(f)) for f in cnt], loc="lower center", ncol=6, frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    fig.text(0.5, 0.995, "Black line: binned median (deciles); grey band: interquartile range; colour: chemical family",
             ha="center", va="top", fontsize=8)
    fig.tight_layout(rect=(0, 0.04, 1, 0.98))
    save(fig, "Fig4_key_relationships")
    out = pd.DataFrame(rows, columns=["target", "feature", "rho_raw", "rho_partial", "ci_lo", "ci_hi", "q"])
    out.to_csv(tab / "table_key_pairs.csv", index=False)
    return out


def step4_group_theory(R, save, tab):
    need = {"n_raman_active_modes", "n_peak", "eta"}
    if not need <= set(R.columns):
        print("  [skip] Fig5: needs n_raman_active_modes")
        return None
    gt = R.dropna(subset=["n_raman_active_modes", "n_peak"])
    gt = gt[gt["n_raman_active_modes"] > 0]
    if len(gt) <= 10:
        print(f"  [skip] Fig5: only {len(gt)} usable rows")
        return None
    cf = fam_colors(R)
    rho = stats.spearmanr(gt["n_raman_active_modes"], gt["n_peak"])[0]
    # Retain the exact analysis rows for the publication-core synthesis figure.
    gt[[c for c in ["mineral", "formula", "family", "n_raman_active_modes", "n_peak", "eta"] if c in gt.columns]] \
        .to_csv(tab / "table_group_theory_ceiling.csv", index=False)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6), gridspec_kw={"width_ratios": [1, 1.2]})
    ax1.scatter(gt["n_raman_active_modes"], gt["n_peak"], s=8, alpha=0.55, lw=0, c=[cf[f] for f in gt["family"]])
    lim = max(gt["n_raman_active_modes"].max(), gt["n_peak"].max())
    ax1.plot([0, lim], [0, lim], "k--", lw=1, label="1:1 (every active mode resolved)")
    ax1.set_xlabel("Group-theoretic Raman-active frequencies, $N_{Raman}$")
    ax1.set_ylabel("Resolved peaks in 100-1300 cm$^{-1}$, $N_{peak}$")
    ax1.text(0.03, 0.97, f"Spearman $\\rho$ = {rho:.2f}\nn = {len(gt)}", transform=ax1.transAxes, va="top", fontsize=8)
    ax1.legend(loc="lower right", frameon=False)
    panel(ax1, "a", dx=-0.1)
    eta = gt[["family", "eta"]].dropna()
    order = eta.groupby("family")["eta"].median().sort_values().index.tolist()
    data = [eta.loc[eta.family == f, "eta"].values for f in order]
    bp = ax2.boxplot(data, patch_artist=True, showfliers=False, widths=0.6, medianprops=dict(color="k", lw=1.2))
    for patch, f in zip(bp["boxes"], order):
        patch.set_facecolor(cf[f])
        patch.set_alpha(0.85)
    ax2.axhline(1, color="k", ls="--", lw=0.8)
    ax2.set_xticks(range(1, len(order) + 1), [f"{fam_name(f)}\n(n={len(d)})" for f, d in zip(order, data)],
                   rotation=40, ha="right")
    ax2.set_ylabel("Mapping efficiency $\\eta = N_{peak}/N_{Raman}$")
    ax2.set_title("$\\eta$ < 1 is expected: $N_{peak}$ counts only 100-1300 cm$^{-1}$", loc="left", fontsize=8.5)
    panel(ax2, "b", dx=-0.08)
    fig.tight_layout()
    save(fig, "Fig5_group_theory_mapping")
    eta_tab = eta.groupby("family")["eta"].agg(["count", "median", lambda s: s.quantile(.25), lambda s: s.quantile(.75)])
    eta_tab.columns = ["n", "median", "q25", "q75"]
    eta_tab.to_csv(tab / "table_eta_by_family.csv")
    return rho


# =============================================================================
# 6. S5 回归 (M1/M2/M3 + 群论子样本)
# =============================================================================
def fit_ols(D, y, preds, cluster_col, fe_family=False, fe_min_n=15):
    preds = [p for p in preds if p in D.columns]
    sub = D.dropna(subset=[y] + preds)
    if len(sub) < len(preds) + 20:
        return None
    X = sub[preds].astype(float)
    keep = [c for c in preds if X[c].std(ddof=0) > 1e-12]
    if not keep:
        return None
    Xz = (X[keep] - X[keep].mean()) / X[keep].std(ddof=0)
    yv = sub[y].astype(float)
    ysd = yv.std(ddof=0)
    yz = ((yv - yv.mean()) / ysd).to_numpy()
    blocks, names = [np.ones((len(sub), 1)), Xz.to_numpy()], ["const"] + keep
    if fe_family:
        fam = sub["family"].where(sub["family"].map(sub["family"].value_counts()) >= fe_min_n, "small_families_pooled")
        dm = pd.get_dummies(fam, drop_first=True).astype(float)
        blocks.append(dm.to_numpy())
        names += [f"FE:{c}" for c in dm.columns]
    Xc = np.hstack(blocks)
    groups = sub[cluster_col].to_numpy() if cluster_col else None
    res = ols_robust(Xc, yz, names, groups)
    res.update(vif=dict(zip(keep, vif(Xz.to_numpy()) if len(keep) > 1 else [1.0])), y_sd=ysd, keep=keep,
               dropped=[c for c in preds if c not in keep])
    return res


def step5_regression(D, tg, cluster_col, save, tab, rep_dir):
    tiers = [("M1 core", CORE_OLS, False), ("M2 extended", EXT_OLS, False), ("M3 extended + family FE", EXT_OLS, True)]
    rows = []
    gt_sub = D.dropna(subset=["n_raman_active_modes"]) if "n_raman_active_modes" in D else D.iloc[:0]
    jobs = [(name, y, preds, fe, D) for y in tg for name, preds, fe in tiers]
    if len(gt_sub) >= 40:
        jobs += [("GT: N_Raman + core", y, ["n_raman_active_modes"] + CORE_OLS[1:], False, gt_sub) for y in GT_RESPONSES if y in tg]
    for name, y, preds, fe, frame in jobs:
        res = fit_ols(frame, y, preds, cluster_col, fe)
        if res is None:
            continue
        for i, t in enumerate(res["names"]):
            if t == "const" or t.startswith("FE:"):
                continue
            rows.append(dict(tier=name, response=y, term=t, beta_std=res["beta"][i], se=res["se"][i],
                             ci_lo=res["lo"][i], ci_hi=res["hi"][i], p=res["p"][i], vif=res["vif"].get(t, np.nan),
                             n=res["n"], r2=res["r2"], adj_r2=res["adj_r2"], se_kind=res["se_kind"],
                             n_clusters=res["n_clusters"], y_sd=res["y_sd"]))
        with open(rep_dir / f"ols_{y}__{name.split(' ')[0].strip(':')}{'_FE' if fe else ''}.txt", "w", encoding="utf-8") as fh:
            fh.write(f"{name}: z({tlabel(y)}) ~ {' + '.join(res['names'][1:])}\n")
            fh.write(f"n={res['n']}, R2={res['r2']:.3f}, adj R2={res['adj_r2']:.3f}, SE={res['se_kind']}"
                     f"{' (G=%d clusters)' % res['n_clusters'] if res['se_kind'] == 'cluster' else ''}\n")
            fh.write("predictors z-scored, response z-scored: coefficients are standardised betas\n\n")
            fh.write(pd.DataFrame({"beta": res["beta"], "se": res["se"], "p": res["p"], "ci_lo": res["lo"], "ci_hi": res["hi"]},
                                  index=res["names"]).round(4).to_string())
            fh.write("\n\nVIF (>5 worrying, >10 severe):\n" + pd.Series(res["vif"]).round(2).to_string() + "\n")
    S = pd.DataFrame(rows)
    if S.empty:
        return S
    S["q"] = S.groupby("tier")["p"].transform(lambda p: bh_fdr(p.values))     # 每个模型层内跨全部系数 BH
    S.to_csv(tab / "table_regression_summary.csv", index=False)

    # Fig6 森林图: M2 vs M3
    resp = [y for y in tg if y in S.response.unique()]
    if resp:
        ncol = 5 if len(resp) > 6 else 3
        nrow = int(np.ceil(len(resp) / ncol))
        fig, axs = plt.subplots(nrow, ncol, figsize=(2.6 * ncol + 1.6, 2.9 * nrow), sharey=True, squeeze=False)
        ypos = {p: i for i, p in enumerate(EXT_OLS[::-1])}
        for ax, y in zip(axs.ravel(), resp):
            for tier, col, off, mk in (("M2 extended", "#0072B2", 0.13, "o"), ("M3 extended + family FE", "#D55E00", -0.13, "D")):
                d = S[(S.response == y) & (S.tier == tier)]
                for _, r in d.iterrows():
                    if r.term not in ypos:
                        continue
                    yy = ypos[r.term] + off
                    filled = r.q < FDR_ALPHA
                    ax.plot([r.ci_lo, r.ci_hi], [yy, yy], color=col, lw=1.1)
                    ax.plot(r.beta_std, yy, mk, ms=4.2, color=col, mfc=col if filled else "white", mew=1)
            ax.axvline(0, color="k", lw=0.7)
            ax.set_title(tlabel(y), fontsize=8.5, loc="left")
            ax.set_xlabel("standardised $\\beta$ (95% CI)")
            ax.set_yticks(range(len(ypos)), [flabel(p) for p in ypos])
            ax.grid(axis="x", color="#eee")
        for ax in axs.ravel()[len(resp):]:
            ax.axis("off")
        fig.legend(handles=[Line2D([], [], marker="o", color="#0072B2", label="M2 extended"),
                            Line2D([], [], marker="D", color="#D55E00", label="M3 extended + family fixed effects"),
                            Line2D([], [], marker="o", color="grey", mfc="white", ls="", label="hollow: BH q >= 0.05")],
                   loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.02))
        fig.tight_layout(rect=(0, 0.04, 1, 1))
        save(fig, "Fig6_regression_forest")
    return S


# =============================================================================
# 7. S6-S7 预测力 / 消融 / LOFO / 置换重要性
# =============================================================================
def design_sets(D, feats):
    sets = {}
    if "raman_snr" in D:
        meas = D[["raman_snr"]].copy()
        if "wl_group" in D:
            meas = pd.concat([meas, pd.get_dummies(D["wl_group"], prefix="wl").astype(float)], axis=1)
        sets["Measurement (SNR, wavelength)"] = meas
    cat = [c for c in ["family", "crystal_system"] if c in D]
    sets["Categorical (family, crystal system)"] = pd.get_dummies(D[cat]).astype(float)
    for g in GROUPS:
        cols = [c for c in feats if fgroup(c) == g]
        if cols:
            sets[GROUP_LABEL[g]] = D[cols]
    sets["All structure/chemistry"] = D[feats]
    return sets


def step6_predictive(D, feats, tg, gcol, args, save, tab):
    sets = design_sets(D, feats)
    ablate = {} if args.fast else {f"All minus {GROUP_LABEL[g]}": D[[c for c in feats if fgroup(c) != g]]
                                   for g in GROUPS if any(fgroup(c) == g for c in feats)}
    all_sets = {**sets, **ablate}
    rows, oof_all, boot_store = [], {}, {}
    t0 = time.time()
    for y in tg:
        m = D[y].notna().to_numpy()
        yy = D.loc[m, y].to_numpy(float)
        grp = D.loc[m, gcol].to_numpy() if gcol else np.arange(m.sum())
        boots = make_boot_idx(grp, args.n_boot, args.seed + 7)
        fold_sets = [group_folds(grp, args.cv_folds, args.seed + r) for r in range(args.repeats)]
        for name, X in all_sets.items():
            Xm = X.loc[m].to_numpy(float)
            preds = np.full((args.repeats, len(yy)), np.nan)
            for r, folds in enumerate(fold_sets):
                for tr, te in folds:
                    preds[r, te] = make_hgb(args.seed, args.fast).fit(Xm[tr], yy[tr]).predict(Xm[te])
            per = [r2_pooled(yy, p) for p in preds]
            bd = boot_r2(yy, preds, boots)
            boot_store[(y, name)] = bd
            lo, hi = np.nanpercentile(bd, [2.5, 97.5])
            rows.append(dict(target=y, feature_set=name, scheme=f"grouped {args.cv_folds}-fold CV x{args.repeats}",
                             r2=np.mean(per), r2_sd_repeats=np.std(per), ci_lo=lo, ci_hi=hi, n=int(m.sum())))
            if name == "All structure/chemistry":
                oof_all[y] = (m, preds)
        print(f"   [CV] {tlabel(y):<28s} done ({time.time() - t0:.0f}s)")
    Rt = pd.DataFrame(rows)

    # 配对差值: All - Categorical / All - Measurement / 消融 dR2
    drows = []
    for y in tg:
        base = boot_store[(y, "All structure/chemistry")]
        ref = Rt[(Rt.target == y) & (Rt.feature_set == "All structure/chemistry")].r2.values[0]

        def add(label, other):
            if (y, other) not in boot_store:
                return
            o = Rt[(Rt.target == y) & (Rt.feature_set == other)].r2.values[0]
            d = base - boot_store[(y, other)]
            drows.append(dict(target=y, contrast=label, delta_r2=ref - o, ci_lo=np.nanpercentile(d, 2.5),
                              ci_hi=np.nanpercentile(d, 97.5)))
        add("All - Categorical", "Categorical (family, crystal system)")
        add("All - Measurement", "Measurement (SNR, wavelength)")
        for g in GROUPS:
            add(f"drop {g}", f"All minus {GROUP_LABEL[g]}")
    Dt = pd.DataFrame(drows)
    Rt.to_csv(tab / "table_cv_performance.csv", index=False)
    Dt.to_csv(tab / "table_cv_paired_differences.csv", index=False)
    return Rt, Dt, oof_all


def step6b_lofo(D, feats, tg, gcol, args, oof_all, tab):
    if args.fast:
        return pd.DataFrame(), pd.DataFrame()
    fsets = {"core": [c for c in CORE_OLS if c in D], "all": feats}
    fam_n = D["family"].value_counts()
    elig = fam_n[fam_n >= args.min_family_size].index.tolist()
    if len(elig) < 2:
        print("  [skip] LOFO: fewer than 2 eligible families")
        return pd.DataFrame(), pd.DataFrame()
    rows, gaps = [], []
    for y in tg:
        m = (D[y].notna() & D["family"].notna()).to_numpy()
        sub = D[m].reset_index(drop=True)
        yy = sub[y].to_numpy(float)
        for fs, cols in fsets.items():
            X = sub[cols].to_numpy(float)
            pred = np.full(len(sub), np.nan)
            base = np.full(len(sub), np.nan)
            for fam in elig:
                te = (sub["family"] == fam).to_numpy()
                tr = ~te
                if te.sum() < 5 or tr.sum() < 30:
                    continue
                pred[te] = make_hgb(args.seed).fit(X[tr], yy[tr]).predict(X[te])
                base[te] = yy[tr].mean()
                sse_b = ((yy[te] - base[te]) ** 2).sum()
                rows.append(dict(feature_set=fs, target=y, held_out_family=fam, n_train=int(tr.sum()), n_test=int(te.sum()),
                                 skill_vs_train_mean=1 - ((yy[te] - pred[te]) ** 2).sum() / sse_b if sse_b > 0 else np.nan,
                                 r2_test_mean=r2_pooled(yy[te], pred[te]),
                                 spearman=stats.spearmanr(yy[te], pred[te])[0] if np.std(pred[te]) > 1e-12 else np.nan))
            ok = np.isfinite(pred)
            sse_b = ((yy[ok] - base[ok]) ** 2).sum()
            pooled = dict(feature_set=fs, target=y, lofo_skill_pooled=1 - ((yy[ok] - pred[ok]) ** 2).sum() / sse_b,
                          lofo_r2_pooled=r2_pooled(yy[ok], pred[ok]), n=int(ok.sum()))
            if fs == "all" and y in oof_all:                     # 同一批(合格族)样本上的随机分组 CV R2
                mm, pr = oof_all[y]
                fam_full = D.loc[mm, "family"].to_numpy()
                yfull = D.loc[mm, y].to_numpy(float)
                e = np.isin(fam_full, elig)
                pooled["grouped_cv_r2_same_rows"] = float(np.mean([r2_pooled(yfull[e], p[e]) for p in pr]))
            gaps.append(pooled)
    L, G = pd.DataFrame(rows), pd.DataFrame(gaps)
    L.to_csv(tab / "table_lofo_by_family.csv", index=False)
    G.to_csv(tab / "table_lofo_pooled_and_gap.csv", index=False)
    return L, G






def step7_importance(D, feats, tg, gcol, args, save, tab, Dt):
    rng = np.random.default_rng(args.seed)
    gidx = {g: [i for i, c in enumerate(feats) if fgroup(c) == g] for g in GROUPS}
    F, G = np.zeros((len(feats), len(tg))), np.zeros((len(GROUPS), len(tg)))
    for j, y in enumerate(tg):
        m = D[y].notna().to_numpy()
        X, yy = D.loc[m, feats].to_numpy(float), D.loc[m, y].to_numpy(float)
        grp = D.loc[m, gcol].to_numpy() if gcol else np.arange(m.sum())
        folds = group_folds(grp, args.cv_folds, args.seed)
        fa, ga = np.zeros(len(feats)), np.zeros(len(GROUPS))
        for tr, te in folds:
            mod = make_hgb(args.seed, args.fast).fit(X[tr], yy[tr])
            Xt, yt = X[te], yy[te]
            v = np.var(yt)
            if v <= 0:
                continue
            base = np.mean((yt - mod.predict(Xt)) ** 2)
            for i in range(len(feats)):
                d = []
                for _ in range(args.n_perm):
                    Xp = Xt.copy()
                    Xp[:, i] = Xp[rng.permutation(len(yt)), i]
                    d.append((np.mean((yt - mod.predict(Xp)) ** 2) - base) / v)
                fa[i] += np.mean(d) / len(folds)
            for k, g in enumerate(GROUPS):
                if not gidx[g]:
                    continue
                d = []
                for _ in range(args.n_perm):
                    Xp = Xt.copy()
                    perm = rng.permutation(len(yt))
                    Xp[:, gidx[g]] = Xt[perm][:, gidx[g]]          # 组内整行联合置换, 保留组内相关
                    d.append((np.mean((yt - mod.predict(Xp)) ** 2) - base) / v)
                ga[k] += np.mean(d) / len(folds)
        F[:, j], G[:, j] = fa, ga
    Fdf = pd.DataFrame(F, index=feats, columns=tg)
    Gdf = pd.DataFrame(G, index=[GROUP_LABEL[g] for g in GROUPS], columns=tg)
    Fdf.to_csv(tab / "table_perm_importance_feature.csv")
    Gdf.to_csv(tab / "table_perm_importance_group.csv")

    order, cols, seps = feat_rows(feats)
    Fo = Fdf.loc[order].clip(lower=0).to_numpy()
    fig = plt.figure(figsize=(12.5, 8.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.3, 1], height_ratios=[1, 1], wspace=0.75, hspace=0.55)
    ax = fig.add_subplot(gs[:, 0])
    im = draw_heat(ax, Fo, None, [tlabel(t) for t in tg], [flabel(c) for c in order], cols, vmin=0,
                   vmax=max(0.05, np.nanpercentile(Fo, 98)), cmap="YlOrRd", fmt="{:.3f}", annot_min=0.008, sep_rows=seps,
                   title="Single-feature permutation importance ($\\Delta R^2$, held-out folds)")
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02).set_label("$\\Delta R^2$")
    panel(ax, "a", dx=-0.6)
    Gv = Gdf.clip(lower=0).to_numpy()
    ax2 = fig.add_subplot(gs[0, 1])
    im2 = draw_heat(ax2, Gv, None, [tlabel(t) for t in tg], list(Gdf.index), [GROUP_COLOR[g] for g in GROUPS], vmin=0,
                    vmax=max(0.1, Gv.max()), cmap="YlOrRd", fmt="{:.2f}", title="Joint group permutation ($\\Delta R^2$)")
    fig.colorbar(im2, ax=ax2, fraction=0.05, pad=0.02)
    panel(ax2, "b", dx=-0.9)
    if Dt is not None and len(Dt):
        A = np.full((len(GROUPS), len(tg)), np.nan)
        Sg = np.full_like(A, np.nan)
        for i, g in enumerate(GROUPS):
            for j, y in enumerate(tg):
                r = Dt[(Dt.target == y) & (Dt.contrast == f"drop {g}")]
                if len(r):
                    A[i, j] = r.delta_r2.values[0]
                    Sg[i, j] = 0.0 if r.ci_lo.values[0] > 0 else 1.0    # CI 不含 0 -> 显著
        ax3 = fig.add_subplot(gs[1, 1])
        im3 = draw_heat(ax3, np.clip(A, 0, None), np.where(Sg == 0, 0.04, 1.0), [tlabel(t) for t in tg],
                        [GROUP_LABEL[g] for g in GROUPS], [GROUP_COLOR[g] for g in GROUPS], vmin=0,
                        vmax=max(0.1, np.nanmax(A)), cmap="YlGnBu", fmt="{:.2f}",
                        title="Drop-group ablation, $\\Delta R^2$ (* = 95% CI excludes 0)")
        fig.colorbar(im3, ax=ax3, fraction=0.05, pad=0.02)
        panel(ax3, "c", dx=-0.9)
    save(fig, "Fig8_importance_and_ablation")
    return Fdf, Gdf


# =============================================================================
# 8. S8 Lasso (嵌套 CV + 稳定性选择) / S9 PLS
# =============================================================================
def lasso_1se(Xs, y, groups, k, seed, n_alpha=60):
    yc = y - y.mean()
    amax = np.max(np.abs(Xs.T @ yc)) / len(y)
    if not np.isfinite(amax) or amax <= 0:
        return 1.0
    alphas = np.geomspace(amax, amax * 1e-3, n_alpha)
    folds = group_folds(groups, k, seed)
    mse = np.zeros((len(folds), n_alpha))
    for i, (tr, va) in enumerate(folds):
        xm, ym = Xs[tr].mean(0), y[tr].mean()
        coefs = lasso_path(Xs[tr] - xm, y[tr] - ym, alphas=alphas, max_iter=20000)[1]
        mse[i] = ((y[va][:, None] - ((Xs[va] - xm) @ coefs + ym)) ** 2).mean(0)
    mu, se = mse.mean(0), mse.std(0, ddof=1) / np.sqrt(len(folds))
    i0 = int(mu.argmin())
    return float(alphas[np.where(mu <= mu[i0] + se[i0])[0].min()])          # 最大的 alpha (最稀疏) 落在 1-SE 内


def lasso_nested_oof(X: pd.DataFrame, y, groups, k, seed):
    pred = np.full(len(y), np.nan)
    for tr, te in group_folds(groups, k, seed):
        pr = Prep().fit(X.iloc[tr])
        Xtr, Xte = pr.transform(X.iloc[tr]), pr.transform(X.iloc[te])
        mu, sd = y[tr].mean(), y[tr].std() or 1.0
        a = lasso_1se(Xtr, (y[tr] - mu) / sd, groups[tr], k, seed + 1)
        pred[te] = Lasso(alpha=a, max_iter=20000).fit(Xtr, (y[tr] - mu) / sd).predict(Xte) * sd + mu
    return pred


def step8_lasso(D, feats, tg, gcol, args, save, tab):
    def lab(n):
        return "Missingness: " + flabel(n[3:]) if n.startswith("NA:") else flabel(n)
    out, coef_rows, freq_rows = {}, {}, {}
    names_ref = None
    for y in tg:
        m = D[y].notna().to_numpy()
        X = D.loc[m, feats].reset_index(drop=True)
        yy = D.loc[m, y].to_numpy(float)
        grp = D.loc[m, gcol].to_numpy() if gcol else np.arange(m.sum())
        pr = Prep().fit(X)
        Xs = pr.transform(X)
        mu, sd = yy.mean(), yy.std()
        ys = (yy - mu) / sd
        alpha = lasso_1se(Xs, ys, grp, args.cv_folds, args.seed)
        coef = Lasso(alpha=alpha, max_iter=20000).fit(Xs, ys).coef_
        # 嵌套 CV R2 (alpha 在每个外折训练集内选择)
        pred = lasso_nested_oof(X, yy, grp, args.cv_folds, args.seed + 3)
        boots = make_boot_idx(grp, args.n_boot, args.seed + 11)
        bd = boot_r2(yy, pred[None, :], boots)
        # bootstrap 稳定性选择
        idx_by_g = [np.asarray(v) for v in pd.Series(np.arange(len(grp))).groupby(grp).apply(list)]
        rng = np.random.default_rng(args.seed + 5)
        cnt = np.zeros(Xs.shape[1])
        for _ in range(args.n_stab):
            idx = np.concatenate([idx_by_g[i] for i in rng.integers(0, len(idx_by_g), len(idx_by_g))])
            cnt += Lasso(alpha=alpha, max_iter=5000).fit(Xs[idx], ys[idx]).coef_ != 0
        freq = cnt / args.n_stab
        names = pr.names
        top = [k for k in np.argsort(-np.abs(coef))[:6] if abs(coef[k]) > 1e-8]
        eq = " ".join(f"{coef[k]:+.2f}*[{lab(names[k])}]({freq[k]:.0%})" for k in top)
        out[y] = dict(r2=r2_pooled(yy, pred), ci=tuple(np.nanpercentile(bd, [2.5, 97.5])), nnz=int((np.abs(coef) > 1e-8).sum()),
                      eq=f"z({tlabel(y)}) ~ {eq}", names=names, alpha=alpha)
        coef_rows[y], freq_rows[y] = pd.Series(coef, index=names), pd.Series(freq, index=names)
        names_ref = names if names_ref is None or len(names) > len(names_ref) else names_ref
        print(f"   [Lasso] {tlabel(y):<28s} nested R2={out[y]['r2']:.2f}  nnz={out[y]['nnz']}")
    C = pd.DataFrame(coef_rows).reindex(index=sorted({n for s in coef_rows.values() for n in s.index})).fillna(0.0)
    Fq = pd.DataFrame(freq_rows).reindex(index=C.index).fillna(0.0)
    C.to_csv(tab / "table_lasso_coefficients.csv")
    Fq.to_csv(tab / "table_lasso_selection_frequency.csv")
    pd.DataFrame({t: {"nested_cv_r2": o["r2"], "ci_lo": o["ci"][0], "ci_hi": o["ci"][1], "n_nonzero": o["nnz"],
                      "alpha_1se": o["alpha"]} for t, o in out.items()}).T.to_csv(tab / "table_lasso_performance.csv")

    order, cols, seps = feat_rows(feats)
    extra = [n for n in C.index if n.startswith("NA:")]
    rows = order + extra
    Cm = C.reindex(rows).to_numpy()
    Fm = Fq.reindex(rows).to_numpy()
    Cm_plot = np.where(np.abs(Cm) < 5e-3, np.nan, Cm)
    vm = np.nanmax(np.abs(Cm_plot)) if np.isfinite(Cm_plot).any() else 0.5
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 8.4), gridspec_kw={"wspace": 0.05})
    xl = [f"{tlabel(t)} ({out[t]['r2']:.2f}, {out[t]['nnz']})" for t in tg]
    ycol = cols + ["#999999"] * len(extra)
    sepr = seps + ([len(order)] if extra else [])
    im = draw_heat(axs[0], Cm_plot, None, xl, [lab(r) for r in rows], ycol, vmin=-vm, vmax=vm, fmt="{:+.2f}", sep_rows=sepr,
                   title="Lasso (1-SE) coefficients; x labels: (nested-CV $R^2$, no. terms)")
    fig.colorbar(im, ax=axs[0], fraction=0.03, pad=0.32, location="bottom", shrink=0.6).set_label("standardised coefficient (grey = 0)")
    im2 = draw_heat(axs[1], np.where(Fm > 0, Fm, np.nan), None, xl, [lab(r) for r in rows], ycol, vmin=0, vmax=1, cmap="Greens",
                    fmt="{:.2f}", show_y=False, sep_rows=sepr, annot_min=0.2,
                    title=f"Selection frequency in {args.n_stab} formula-cluster bootstraps")
    fig.colorbar(im2, ax=axs[1], fraction=0.03, pad=0.32, location="bottom", shrink=0.6).set_label("selection frequency")
    panel(axs[0], "a", dx=-0.55)
    panel(axs[1], "b", dx=-0.02)
    group_legend(fig, loc="upper center", bbox_to_anchor=(0.5, 0.985))
    save(fig, "Fig9_lasso_sparse_relations")
    return out, C, Fq


def step9_pls(D, feats, tg, gcol, args, save, tab):
    ytg = [t for t in tg if t not in PLS_EXCLUDE]
    m = D[ytg].notna().all(axis=1).to_numpy()
    X, Y = D.loc[m, feats].reset_index(drop=True), D.loc[m, ytg].to_numpy(float)
    grp = D.loc[m, gcol].to_numpy() if gcol else np.arange(m.sum())
    q2 = []
    for k in range(1, 7):
        press = sst = 0.0
        for tr, te in group_folds(grp, args.cv_folds, args.seed):
            pr = Prep().fit(X.iloc[tr])
            mu, sd = Y[tr].mean(0), Y[tr].std(0)
            sd[sd == 0] = 1
            p = PLSRegression(n_components=k, scale=False).fit(pr.transform(X.iloc[tr]), (Y[tr] - mu) / sd) \
                .predict(pr.transform(X.iloc[te]))
            Yt = (Y[te] - mu) / sd
            press += ((Yt - p) ** 2).sum()
            sst += (Yt ** 2).sum()                                     # 相对训练均值(=0)的总平方和
        q2.append(1 - press / sst)
    pd.DataFrame({"n_components": range(1, 7), "Q2_cumulative": q2}).to_csv(tab / "table_pls_q2.csv", index=False)
    pr = Prep().fit(X)
    Xs = pr.transform(X)
    Yz = (Y - Y.mean(0)) / Y.std(0)
    pls = PLSRegression(n_components=2, scale=False).fit(Xs, Yz)
    T = pls.x_scores_
    cf = fam_colors(D)
    fam = D.loc[m, "family"].to_numpy()
    fig = plt.figure(figsize=(12.5, 5.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1, 0.85], wspace=0.5)
    ax = fig.add_subplot(gs[0])
    ax.scatter(T[:, 0], T[:, 1], s=7, c=[cf[f] for f in fam], alpha=0.55, lw=0)
    ax.legend(handles=[Patch(color=cf[f], label=fam_name(f)) for f in pd.Series(fam).value_counts().index[:8]],
              fontsize=6.5, frameon=False, loc="upper left", ncol=2)
    ax.axhline(0, color="#aaa", lw=0.6)
    ax.axvline(0, color="#aaa", lw=0.6)
    ax.set_xlabel("Structural latent variable LV1")
    ax.set_ylabel("Structural latent variable LV2")
    ax.set_title(f"Scores (cumulative $Q^2$: 1 LV = {q2[0]:.2f}, 2 LV = {q2[1]:.2f})", loc="left")
    panel(ax, "a", dx=-0.1)
    ax2 = fig.add_subplot(gs[1])
    order, cols, _ = feat_rows(feats)
    li = {c: i for i, c in enumerate(pr.names)}
    ld = np.array([pls.x_loadings_[li[c], 0] for c in order])
    ld2 = np.array([pls.x_loadings_[li[c], 1] for c in order])
    o = np.argsort(ld)
    yp = np.arange(len(o))
    ax2.barh(yp - 0.18, ld[o], 0.36, color=[cols[i] for i in o])
    ax2.barh(yp + 0.18, ld2[o], 0.36, color=[cols[i] for i in o], alpha=0.4, hatch="////", edgecolor="white")
    ax2.set_yticks(yp, [flabel(order[i]) for i in o], fontsize=6.6)
    ax2.axvline(0, color="k", lw=0.7)
    ax2.set_title("X loadings (solid: LV1; hatched: LV2)", loc="left")
    panel(ax2, "b", dx=-0.02)
    ax3 = fig.add_subplot(gs[2])
    yl = pls.y_loadings_
    yy = np.arange(len(ytg))
    ax3.barh(yy - 0.2, yl[:, 0], 0.4, color="#333", label="LV1")
    ax3.barh(yy + 0.2, yl[:, 1], 0.4, color="#E69F00", label="LV2")
    ax3.set_yticks(yy, [tlabel(t) for t in ytg])
    ax3.invert_yaxis()
    ax3.axvline(0, color="k", lw=0.7)
    ax3.legend(frameon=False)
    ax3.set_title("Spectral loadings", loc="left")
    panel(ax3, "c", dx=-0.02)
    group_legend(fig, loc="lower center", bbox_to_anchor=(0.5, -0.06))
    save(fig, "Fig10_pls_latent_structure")
    pd.DataFrame({"feature": order, "LV1": [pls.x_loadings_[li[c], 0] for c in order],
                  "LV2": [pls.x_loadings_[li[c], 1] for c in order]}).to_csv(tab / "table_pls_x_loadings.csv", index=False)
    return q2


# =============================================================================
# 9. S10 稳健性
# =============================================================================
def step10_sensitivity(raw_df, args, feats, tg, covars, tab):
    specs = [(args.min_snr, args.min_peaks, "main")] + [(s, args.min_peaks, f"SNR>={s:g}") for s in (10, 40) if s != args.min_snr] \
            + [(args.min_snr, p, f"peaks>={p}") for p in (2, 5) if p != args.min_peaks]
    tabs = {}
    for snr, pk, lab in specs:
        d, _ = apply_qc(raw_df, args.qualities, snr, pk, args.orientation, verbose=False)
        _, D2, _ = build_frames(d)
        tabs[lab] = (corr_table(D2, feats, tg, covars), len(D2))
    main = tabs["main"][0].set_index(["feature", "target"])
    strong = main[(main.q < FDR_ALPHA) & (main.rho.abs() >= 0.2)].index          # 主分析中"强且显著"的配对
    rows = []
    for lab, (T, n) in tabs.items():
        if lab == "main":
            continue
        T = T.set_index(["feature", "target"])
        Ts, Ms = T.reindex(strong), main.loc[strong]
        both = pd.concat([main.rho, T.rho.reindex(main.index)], axis=1, keys=["m", "s"]).dropna()
        rows.append(dict(spec=lab, n=n, n_strong_pairs_main=len(strong),
                         frac_same_sign=float((np.sign(Ts.rho) == np.sign(Ms.rho)).mean()) if len(strong) else np.nan,
                         frac_same_sign_and_q_lt_0p05=float(((np.sign(Ts.rho) == np.sign(Ms.rho)) & (Ts.q < FDR_ALPHA)).mean()) if len(strong) else np.nan,
                         spearman_of_all_rho_vs_main=float(stats.spearmanr(both.m, both.s)[0]) if len(both) > 5 else np.nan))
    S = pd.DataFrame(rows)
    S.to_csv(tab / "table_sensitivity_qc.csv", index=False)
    return S


def md_table(df: pd.DataFrame) -> str:
    if df is None or df.empty:
        return "_(none)_"
    head = "| " + " | ".join(map(str, df.columns)) + " |\n|" + "---|" * len(df.columns) + "\n"
    body = "\n".join("| " + " | ".join("" if (isinstance(v, float) and not np.isfinite(v)) else
                                       (f"{v:.2f}" if isinstance(v, (float, np.floating)) else str(v)) for v in r) + " |"
                     for r in df.itertuples(index=False))
    return head + body



# =============================================================================
# 10b. 投稿版 Figure pipeline: 7 main + 6 supplementary figures
# =============================================================================
# =============================================================================
# 10b. 投稿版 Figure pipeline: 7 main + 6 supplementary  (v2: 版式校验 / 数据驱动文字)
# =============================================================================
SHORT_T = {"n_peak": "ln N_peak", "evenness": "Evenness", "gamma": "ln \u03b3", "domega_median": "ln \u0394\u03c9",
           "r_overlap": "Overlap", "frac_low": "Low frac.", "frac_mid": "Mid frac.", "frac_high": "High frac.",
           "w1_distance": "ln W1", "eta": "\u03b7"}
SHORT_FAM = {"phosphate_arsenate_vanadate": "phosphate/As/V", "sulfide_selenide_telluride": "sulfide/Se/Te",
             "oxide_hydroxide": "oxide/hydroxide", "native_element": "native element"}
SET_SHORT = [("Measurement (SNR, wavelength)", "Measurement"), ("Categorical (family, crystal system)", "Family + system"),
             ("Size / complexity", "Size"), ("Symmetry / group theory", "Symmetry"), ("Bond geometry", "Bond geometry"),
             ("Chemistry", "Comp.-derived chemistry"), ("Local environment (SOAP)", "SOAP"), ("All structure/chemistry", "All descriptors")]

# 先验假设驱动的 headline 配对: 必须在看结果之前固定. 若是看过结果后才选的, 论文中须声明为探索性.
# eta = N_peak/N_Raman 与 N_Raman 相关的描述符(raman_mode_fraction, degeneracy_fraction, ...)存在定义性耦合, 故不作 headline.
HEADLINE_PAIRS = [("mean_reduced_mass", "frac_low"), ("mean_reduced_mass", "n_peak"), ("mean_reduced_mass", "frac_high"),
                  ("mean_delta_chi_bond", "w1_distance"), ("mixing_entropy", "gamma"), ("sigma_chi", "frac_high"),
                  ("n_anion_group_types", "domega_median"), ("n_anion_group_types", "n_peak")]
HEADLINE_MULTI = [("gamma", "mean_reduced_mass"), ("gamma", "sigma_chi"), ("gamma", "mixing_entropy"),
                  ("n_peak", "mean_reduced_mass"), ("frac_low", "mean_reduced_mass"),
                  ("w1_distance", "mean_delta_chi_bond")]


def tshort(c):
    return SHORT_T.get(c, tlabel(c))


def fshort(f):
    return SHORT_FAM.get(f, fam_name(f))


def _finalize_figure(fig, layout=True):
    """落实面板字母, 隐藏总标题, 按轴宽折行标题, 启用 constrained layout."""
    import textwrap
    for ax in fig.axes:
        _apply_panel(ax)
    lettered = [ax for ax in fig.axes if getattr(ax, "_panel_done", False)]
    if len(lettered) == 1:                             # 单面板图不需要 (a)
        ax = lettered[0]
        t, pre = ax.get_title(loc="left"), f"({ax._panel_letter}) "
        if t.startswith(pre):
            ax.set_title(t[len(pre):], loc="left")
        art = getattr(ax, "_panel_artist", None)
        if art is not None:
            art.remove()
    st = getattr(fig, "_suptitle", None)



    if st is not None and not FIG_TITLES:
        st.remove()
        fig._suptitle = None
    if layout:
        try:
            fig.set_layout_engine("constrained")
            fig.canvas.draw()                              # 先算出真实轴宽, 再折行标题
        except Exception:
            pass
    w_in = fig.get_size_inches()[0]
    for ax in fig.axes:
        t = ax.get_title(loc="left")
        if t and "\n" not in t:
            ax.set_title(textwrap.fill(t, max(16, int(ax.get_position().width * w_in * 12.5))), loc="left")


def figure_qa(fig, name, overlap_frac=0.15):
    """渲染后检查: (1) MathText 逐行可解析 (2) 任意两个文字对象的包围盒不重叠. 返回问题列表."""
    import itertools
    from matplotlib.mathtext import MathTextParser
    from matplotlib.text import Text
    issues = []
    parser = MathTextParser("agg")
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    hidden = set()                                     # 超出坐标范围、实际不会绘制的刻度标签
    for ax in fig.axes:
        if not ax.axison:                              # axis("off"): 刻度标签不绘制
            hidden |= set(ax.get_xticklabels()) | set(ax.get_yticklabels())
            continue
        try:
            (xa, xb), (ya, yb) = sorted(ax.get_xlim()), sorted(ax.get_ylim())
        except Exception:
            continue
        tx, ty = 1e-6 * (xb - xa), 1e-6 * (yb - ya)
        hidden |= {l for l in ax.get_xticklabels() if not (xa - tx <= l.get_position()[0] <= xb + tx)}
        hidden |= {l for l in ax.get_yticklabels() if not (ya - ty <= l.get_position()[1] <= yb + ty)}
    items = []
    for t in fig.findobj(Text):
        s = t.get_text()
        if not t.get_visible() or not s.strip() or t in hidden:
            continue
        for line in s.split("\n"):                      # 逐行解析: 多行文本整体送入 MathText 会误报 \n
            if line.count("$") >= 2:
                if "\\\\" in line:
                    issues.append(f"doubled LaTeX backslash: {line[:60]!r}")
                try:
                    parser.parse(line, dpi=fig.dpi)
                except Exception as exc:
                    issues.append(f"MathText failure {line[:60]!r}: {exc}")
        try:
            bb = t.get_window_extent(rend)
        except Exception:
            continue
        if bb.width >= 1 and bb.height >= 1:
            items.append((t, bb))
    n_ov = 0
    for (a, ba), (b, bb) in itertools.combinations(items, 2):
        if a.get_rotation() % 90 and b.get_rotation() % 90:      # 平行斜排标签: 轴对齐包围盒相交不代表字形相交
            continue
        w = min(ba.x1, bb.x1) - max(ba.x0, bb.x0)
        h = min(ba.y1, bb.y1) - max(ba.y0, bb.y0)
        if w > 0 and h > 0 and w * h > overlap_frac * min(ba.width * ba.height, bb.width * bb.height):
            n_ov += 1
            if n_ov <= 4:
                issues.append(f"text overlap: {a.get_text()[:32]!r} x {b.get_text()[:32]!r}".replace("\n", " "))
    if n_ov > 4:
        issues.append(f"... {n_ov} overlapping text pairs in total")
    QA_LOG.append((name, n_ov, issues))
    return issues


def _save_fig(fig, out, sub, name, layout=True):
    _finalize_figure(fig, layout)
    issues = figure_qa(fig, name)
    if issues:
        print(f"  [figure QA] {name}: {len(issues)} issue(s); first: {issues[0]}")
    Saver(out / "figures" / sub)(fig, name)


def _save_pubfig(fig, out, name, layout=True):
    _save_fig(fig, out, "main", name, layout)


def _save_suppfig(fig, out, name, layout=True):
    _save_fig(fig, out, "supplementary", name, layout)


def _fmt_ci(r, digits=2):
    return f"{r.rho:+.{digits}f} [{r.ci_lo:+.{digits}f}, {r.ci_hi:+.{digits}f}]"


def _headline_pairs(par, n=8):
    """先验固定的 headline 配对; 不足 n 个时按 |rho| 补足 (排除 eta 与定义性配对)."""
    rows = []
    for x, y in HEADLINE_PAIRS:
        z = par[(par.feature == x) & (par.target == y)]
        if len(z) and np.isfinite(z.iloc[0].rho):
            rows.append(z.iloc[0])
    if len(rows) < n:
        used = {(r.feature, r.target) for r in rows}
        extra = par.dropna(subset=["rho"]).assign(a=lambda d: d.rho.abs()).sort_values("a", ascending=False)
        for _, r in extra.iterrows():
            if r.target in ML_EXCLUDE or (r.feature, r.target) in used or (r.feature, r.target) in TRIVIAL_PAIRS:
                continue
            rows.append(r)
            if len(rows) >= n:
                break
    return pd.DataFrame(rows).head(n)


def _evidence_sentences(par, Rt, Dt, Gdf, rho_gt, sens):
    """全部由本次运行的数值生成; 不含预先写好的机制结论."""
    L = []
    if rho_gt is not None:
        s = "weak" if rho_gt < 0.4 else "moderate" if rho_gt < 0.7 else "strong"
        L.append(f"Group-theoretical N_Raman and observed N_peak show a {s} monotonic association (Spearman rho = {rho_gt:.2f}).")
    if Dt is not None and len(Dt):
        d = Dt[Dt.contrast == "All - Categorical"]
        yes = d[d.ci_lo > 0]
        no = d[d.ci_lo <= 0]
        if len(yes):
            L.append("Structural and composition-derived descriptors add out-of-fold information beyond family + crystal system (paired dR2 95% CI > 0) for: "
                     + "; ".join(f"{tlabel(r.target)} ({r.delta_r2:+.2f} [{r.ci_lo:+.2f}, {r.ci_hi:+.2f}])" for r in yes.itertuples()) + ".")
        if len(no):
            L.append("No detectable gain over family + crystal system (CI includes or is below 0) for: "
                     + ", ".join(tlabel(t) for t in no.target) + ".")
        m = Dt[Dt.contrast == "All - Measurement"]
        if len(m) and Rt is not None and len(Rt):
            meas = Rt[Rt.feature_set == "Measurement (SNR, wavelength)"].set_index("target")
            # 仅当 Measurement 自身有正的样本外解释力 (CI 下界 > 0) 且 All 明显更差时, 才称"仪器主导"
            worse = m[(m.ci_hi < 0) & m.target.map(lambda t: t in meas.index and meas.loc[t, "ci_lo"] > 0)]
            if len(worse):
                L.append("Measurement conditions alone out-predict all structural and composition-derived descriptors for: "
                         + ", ".join(tlabel(t) for t in worse.target) + " (treat these endpoints as instrument-dominated).")
    if Rt is not None and len(Rt):
        a = Rt[Rt.feature_set == "All structure/chemistry"]
        pos = a[a.ci_lo > 0]
        L.append("Endpoints with out-of-fold R2 whose 95% CI excludes 0: "
                 + (", ".join(f"{tlabel(r.target)} (R2={r.r2:.2f})" for r in pos.itertuples()) if len(pos) else "none") + ".")
    if Gdf is not None and len(Gdf):
        parts = []
        for t in Gdf.columns:
            col = Gdf[t]
            if col.max() >= 0.02:
                parts.append(f"{tlabel(t)}: {col.idxmax()} (dR2={col.max():.2f})")
        L.append("Largest joint-permutation group per endpoint (dR2 >= 0.02): " + ("; ".join(parts) if parts else "none reached 0.02") + ".")
    if sens is not None and len(sens):
        L.append(f"QC sensitivity across {len(sens)} alternative settings: sign agreement of main strong partial correlations "
                 f"{sens.frac_same_sign.min():.2f}-{sens.frac_same_sign.max():.2f}; sign agreement with q<0.05 retained "
                 f"{sens.frac_same_sign_and_q_lt_0p05.min():.2f}-{sens.frac_same_sign_and_q_lt_0p05.max():.2f}; "
                 f"rank agreement of the whole partial-correlation matrix {sens.spearman_of_all_rho_vs_main.min():.2f}-"
                 f"{sens.spearman_of_all_rho_vs_main.max():.2f}.")
    return L


def write_paper_results(out, args, n, par, S, Rt, Dt, Gdf, rho_gt, sens):
    rows = []
    hp = _headline_pairs(par)
    for _, r in hp.iterrows():
        rows.append(dict(section="Association", structure=flabel(r.feature), spectrum=tlabel(r.target),
                         estimate=r.rho, ci_lo=r.ci_lo, ci_hi=r.ci_hi, q=r.q, n=int(r.n)))
    if S is not None and len(S):
        for response, term in HEADLINE_MULTI:
            for tier in ("M3 extended + family FE", "M2 extended"):
                q = S[(S.tier == tier) & (S.response == response) & (S.term == term)]
                if len(q):
                    r = q.iloc[0]
                    rows.append(dict(section=f"Multivariable ({tier})", structure=flabel(term), spectrum=tlabel(response),
                                     estimate=r.beta_std, ci_lo=r.ci_lo, ci_hi=r.ci_hi, q=r.q, n=int(r.n)))
    if Rt is not None and len(Rt):
        for r in Rt[Rt.feature_set == "All structure/chemistry"].itertuples():
            rows.append(dict(section="Grouped-CV R2", structure="All structure/chemistry", spectrum=tlabel(r.target),
                             estimate=r.r2, ci_lo=r.ci_lo, ci_hi=r.ci_hi, q=np.nan, n=int(r.n)))
    if Dt is not None and len(Dt):
        for r in Dt[Dt.contrast.isin(["All - Categorical", "All - Measurement"])].itertuples():
            rows.append(dict(section="Incremental dR2", structure=r.contrast, spectrum=tlabel(r.target),
                             estimate=r.delta_r2, ci_lo=r.ci_lo, ci_hi=r.ci_hi, q=np.nan, n=np.nan))
    pd.DataFrame(rows).to_csv(out / "tables" / "table_paper_key_results.csv", index=False)
    lines = ["# Data-driven results statements\n",
             "Every sentence below is generated from this run's numbers. Mechanistic interpretation is intentionally NOT pre-written: "
             "add it to the manuscript only where the statements support it.\n",
             f"Headline association pairs were fixed in code (HEADLINE_PAIRS) and are hypothesis-driven; if they were chosen after "
             f"inspecting earlier results, report them as exploratory. n = {n}.\n"]
    lines += [f"- {s}" for s in _evidence_sentences(par, Rt, Dt, Gdf, rho_gt, sens)]
    (out / "paper_results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---- captions: purely descriptive (no result claims); results statements live in paper_results.md ----
MAIN_CAPTIONS = {
    "Figure 1": "Study design. Raman spectra from the RRUFF Raman directory are linked to crystallographic and composition information from the RRUFF AMCSD directory; structural and composition-derived descriptors are then constructed for prediction. (a) Analysis framework; (b) quality-control funnel; (c) median standardised spectral profile per chemical family.",
    "Figure 2": "Group-theoretical mode availability versus observed peaks. (a) N_Raman against N_peak with the 1:1 line (all active modes resolved). (b, c) Partial Spearman association of symmetry descriptors with the peak-count residual given ln N_Raman, i.e. resolved peaks beyond what mode availability implies; the ratio eta = N_peak/N_Raman is deliberately not used because it is coupled to N_Raman by construction.",
    "Figure 3": "Bond-geometry descriptors against spectral endpoints. Line: decile-binned median. Annotations: partial Spearman rho (adjusted for ln N_atom and ln SNR), 95% CI and BH-FDR q.",
    "Figure 4": "Composition-derived chemistry descriptors against spectral endpoints (same conventions as Figure 3).",
    "Figure 5": "Out-of-fold R2 (grouped 5-fold CV by chemical formula) for each predictor set and spectral endpoint. Negative values (blue) are worse than predicting the mean. Measurement = SNR + excitation wavelength; Family + system = chemical family and crystal system.",
    "Figure 6": "Joint held-out permutation importance (change in R2 when an entire descriptor group is permuted jointly). Values quantify predictive information, not causal effect.",
    "Figure 7": "Proposed conceptual framework (hypothesis). Boxes are the analysed layers; arrows indicate the hypothesised direction and are tested only through the associations reported in Figures 2-6.",
}
SUPP_CAPTIONS = {
    "Figure S1": "Quality-control funnel and median standardised spectral profile per chemical family.",
    "Figure S2": "Complete partial Spearman correlation matrix (adjusted for ln N_atom, ln SNR). * q<0.05, ** q<0.01, *** q<0.001 (BH-FDR over all cells).",
    "Figure S3": "Kruskal-Wallis effect size eta^2_H of categorical factors on each spectral endpoint.",
    "Figure S4": "Standardised OLS coefficients (95% CI, formula-cluster-robust SE) for the extended model (M2) and the extended model with chemical-family fixed effects (M3). Filled: BH q<0.05.",
    "Figure S5": "Leave-one-family-out skill score relative to the training-set mean for each held-out family.",
    "Figure S6": "(a) Lasso coefficients (top 18 by magnitude; 'Missingness' rows are indicator variables). (b) Cumulative grouped-CV Q2 of PLS. (c) Rank agreement of the partial-correlation matrix under alternative QC thresholds.",
}


# ---- plotting helpers ----
def _scatter_binned(ax, D, par, x, y):
    if x not in D or y not in D:
        ax.axis("off")
        return
    q = D[[x, y]].dropna()
    ax.scatter(q[x], q[y], s=6, alpha=0.30, lw=0, color="0.35")
    if len(q) > 40 and q[x].nunique() > 10:
        bins = pd.qcut(q[x], 10, duplicates="drop")
        med = q.groupby(bins, observed=True).median(numeric_only=True)
        ax.plot(med[x], med[y], "-", color="#D55E00", lw=1.6)
    z = par[(par.feature == x) & (par.target == y)]
    if len(z) and np.isfinite(z.iloc[0].rho):
        z = z.iloc[0]
        txt = f"partial $\\rho$ = {z.rho:+.2f}\n[{z.ci_lo:+.2f}, {z.ci_hi:+.2f}]\nq {fmt_p(z.q)}"
        if z.rho >= 0:
            annotate_stats(ax, txt)
        else:
            annotate_stats(ax, txt, x=0.97, ha="right")
    ax.set_xlabel(flabel(x))
    ax.set_ylabel(tlabel(y))


def _peak_residual(D):
    """ln N_peak 对 ln N_Raman 回归后的残差: 超出'模式可用数'所隐含的可分辨峰数."""
    if not {"n_peak", "n_raman_active_modes"} <= set(D.columns):
        return None
    q = D[["n_peak", "n_raman_active_modes"]].dropna()
    if len(q) < 40:
        return None
    X = np.column_stack([np.ones(len(q)), q["n_raman_active_modes"]])
    b = np.linalg.lstsq(X, q["n_peak"], rcond=None)[0]
    return pd.Series(q["n_peak"].to_numpy() - X @ b, index=q.index)


def _forest_grid(S, resp, figsize=(7.15, 7.4), ncol=3):
    nrow = int(np.ceil(len(resp) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=figsize, sharey=True, squeeze=False)
    order = EXT_OLS[::-1]
    ypos = {p: i for i, p in enumerate(order)}
    for k, (ax, y) in enumerate(zip(axs.ravel(), resp)):
        for tier, col, off, mk in (("M2 extended", "#0072B2", 0.14, "o"), ("M3 extended + family FE", "#D55E00", -0.14, "D")):
            for r in S[(S.response == y) & (S.tier == tier)].itertuples():
                if r.term not in ypos:
                    continue
                yy = ypos[r.term] + off
                ax.plot([r.ci_lo, r.ci_hi], [yy, yy], color=col, lw=0.9)
                ax.plot(r.beta_std, yy, mk, ms=3.6, color=col, mfc=col if r.q < FDR_ALPHA else "white", mew=0.9)
        ax.axvline(0, color="k", lw=0.6)
        ax.set_title(tlabel(y), loc="left")
        ax.set_xlabel("Standardised $\\beta$")
        ax.set_yticks(range(len(order)), [flabel(p) for p in order])
        ax.tick_params(axis="y", labelleft=(k % ncol == 0))
        ax.grid(axis="x", color="0.92")
    for ax in axs.ravel()[len(resp):]:
        ax.axis("off")
    fig.legend(handles=[Line2D([], [], marker="o", color="#0072B2", label="M2 extended"),
                        Line2D([], [], marker="D", color="#D55E00", label="M3 + family fixed effects"),
                        Line2D([], [], marker="o", color="grey", mfc="white", ls="", label="hollow: q >= 0.05")],
               loc="outside lower center", ncol=3, frameon=False)
    return fig


def publication_figure_pipeline(out, args, R, D, funnel, par, cat, S, Rt, Dt, L, G, lasso, C, Fq, q2, Gdf, rho_gt, sens, tg):
    """7 main + 6 supplementary figures. 不改变任何分析结果; 版式经 figure_qa 渲染检查."""
    import textwrap
    out = Path(out)
    for sub in ("main", "supplementary"):
        (out / "figures" / sub).mkdir(parents=True, exist_ok=True)
    W = 7.15                                              # 双栏图宽 (英寸)

    # # ---------------- Figure 1 ----------------
    # fig = plt.figure(figsize=(W, 6.2), layout="constrained")
    # gs = fig.add_gridspec(2, 2, height_ratios=[0.95, 1.05], width_ratios=[0.85, 1.25])
    # ax = fig.add_subplot(gs[0, :]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    # panel(ax, "a", dx=0.0, dy=0.98)
    # boxes = [
    #     (0.04, 0.62, 0.40, 0.26, "RRUFF AMCSD\ncrystallographic + composition records"),
    #     (0.72, 0.62, 0.22, 0.26, "RRUFF Raman\ndirectory"),
    #     (0.01, 0.12, 0.19, 0.26, "Symmetry /\ngroup theory"),
    #     (0.24, 0.12, 0.19, 0.26, "Bond\ngeometry"),
    #     (0.47, 0.12, 0.19, 0.26, "Composition-derived\nchemistry"),
    #     (0.74, 0.12, 0.23, 0.26, "Spectral endpoints\n(QC cohort)")
    # ]
    # for x, y, w, h, t in boxes:
    #     ax.add_patch(Rectangle((x, y), w, h, fill=False, lw=0.9, ec="0.25"))
    #     ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=8)
    # arr = dict(arrowstyle="->", lw=0.9)
    # for a0, b0 in [
    #     ((0.14, 0.62), (0.105, 0.38)),
    #     ((0.29, 0.62), (0.31, 0.38)),
    #     ((0.45, 0.62), (0.56, 0.38)),
    #     ((0.83, 0.62), (0.855, 0.38))
    # ]:
    #     ax.annotate("", xy=b0, xytext=a0, arrowprops=arr)
    # for x0 in (0.20, 0.43, 0.66):                                                        # descriptor layers vs endpoints (tested)
    #     ax.annotate("", xy=(0.74, 0.25), xytext=(x0, 0.25), arrowprops=dict(arrowstyle="->", lw=0.9, ls="--", color="0.45"))
    # ax.text(0.5, 0.03, "dashed arrows: associations tested in this study", ha="center", fontsize=7, color="0.35")
    # ax = fig.add_subplot(gs[1, 0])
    # y = np.arange(len(funnel))[::-1]
    # ax.barh(y, funnel["n"], color="0.55", height=0.62)
    # for yi, nv in zip(y, funnel["n"]):
    #     ax.text(nv + funnel["n"].max() * 0.02, yi, f"{int(nv)}", va="center", fontsize=7)
    # ax.set_yticks(y, [textwrap.fill(s, 20) for s in funnel["step"]])
    # ax.set_xlabel("Number of spectra"); ax.set_xlim(0, funnel["n"].max() * 1.2)
    # ax.set_title("QC funnel", loc="left"); panel(ax, "b")
    # ax = fig.add_subplot(gs[1, 1])
    # z = (D[tg] - D[tg].mean()) / D[tg].std(); z["family"] = D["family"]
    # prof = z.groupby("family").median(numeric_only=True); cnt = D["family"].value_counts(); prof = prof.loc[cnt.index]
    # draw_heat(ax, prof.values, None, [tshort(t) for t in tg], [f"{fshort(f)} ({cnt[f]})" for f in prof.index],
    #           vmin=-1.2, vmax=1.2, fmt="{:.1f}", fs=5.8, title="Median spectral profile (z-score)")
    # panel(ax, "c")
    # _save_pubfig(fig, out, "Figure1_framework_and_cohort")


    # ---------------- Figure 1 ----------------
    fig = plt.figure(figsize=(W, 6.2), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[0.95, 1.05], width_ratios=[0.85, 1.25])
    ax = fig.add_subplot(gs[0, :]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_title("Framework and cohort", loc="left");
    panel(ax, "a", dx=0.0, dy=0.98)

    boxes = [
        (0.03, 0.58, 0.45, 0.30, "RRUFF AMCSD\ncrystallographic + composition records"),
        (0.52, 0.58, 0.45, 0.30, "RRUFF Raman\nspectral library"),
        (0.03, 0.06, 0.55, 0.30, "Symmetry / group theory + bond geometry \n + composition-derived chemistry"),
        (0.68, 0.06, 0.29, 0.30, "Spectral endpoints\n(QC cohort)")
    ]
    for x, y, w, h, t in boxes:
        ax.add_patch(Rectangle((x, y), w, h, fill=False, lw=0.9, ec="0.25"))
        ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=8)

    arr = dict(arrowstyle="->", lw=0.9)
    # AMCSD -> descriptor block
    ax.annotate("", xy=(0.255, 0.36), xytext=(0.255, 0.58), arrowprops=arr)
    # Raman -> spectral endpoints
    ax.annotate("", xy=(0.82, 0.36), xytext=(0.82, 0.58), arrowprops=arr)
    # descriptor block -> spectral endpoints
    ax.annotate("", xy=(0.68, 0.21), xytext=(0.58, 0.21), arrowprops=arr)

    ax = fig.add_subplot(gs[1, 0])
    y = np.arange(len(funnel))[::-1]
    ax.barh(y, funnel["n"], color="0.55", height=0.62)
    for yi, nv in zip(y, funnel["n"]):
        ax.text(nv + funnel["n"].max() * 0.02, yi, f"{int(nv)}", va="center", fontsize=7)
    ax.set_yticks(y, [textwrap.fill(s, 20) for s in funnel["step"]])
    ax.set_xlabel("Number of spectra"); ax.set_xlim(0, funnel["n"].max() * 1.2)
    ax.set_title("QC funnel", loc="left"); panel(ax, "b")
    ax = fig.add_subplot(gs[1, 1])
    z = (D[tg] - D[tg].mean()) / D[tg].std(); z["family"] = D["family"]
    prof = z.groupby("family").median(numeric_only=True); cnt = D["family"].value_counts(); prof = prof.loc[cnt.index]
    draw_heat(ax, prof.values, None, [tshort(t) for t in tg], [f"{fshort(f)} ({cnt[f]})" for f in prof.index],
            vmin=-1.2, vmax=1.2, fmt="{:.1f}", fs=5.8, title="Median spectral profile (z-score)")
    panel(ax, "c")
    _save_pubfig(fig, out, "Figure1_framework_and_cohort")


    # # ---------------- Figure 1 ----------------
    # fig = plt.figure(figsize=(W, 6.2), layout="constrained")
    # gs = fig.add_gridspec(2, 2, height_ratios=[0.95, 1.05], width_ratios=[0.85, 1.25])
    # ax = fig.add_subplot(gs[0, :]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    # ax.set_title("Framework  Overview", loc="left");
    # panel(ax, "a", dx=0.0, dy=0.98)

    # boxes = [
    #     (0.03, 0.60, 0.44, 0.28, "RRUFF AMCSD\ncrystallographic + composition records"),
    #     (0.55, 0.60, 0.42, 0.28, "RRUFF Raman\ndirectory"),
    #     (0.02, 0.08, 0.17, 0.28, "Symmetry /\ngroup theory"),
    #     (0.21, 0.08, 0.17, 0.28, "Bond\ngeometry"),
    #     (0.40, 0.08, 0.17, 0.28, "Composition-derived\nchemistry"),
    #     (0.60, 0.08, 0.38, 0.28, "Spectral endpoints\n(QC cohort)")
    # ]
    # for x, y, w, h, t in boxes:
    #     ax.add_patch(Rectangle((x, y), w, h, fill=False, lw=0.9, ec="0.25"))
    #     ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=8)

    # arr = dict(arrowstyle="->", lw=0.9)
    # for a0, b0 in [
    #     ((0.25, 0.60), (0.105, 0.36)),
    #     ((0.25, 0.60), (0.295, 0.36)),
    #     ((0.25, 0.60), (0.485, 0.36)),
    #     ((0.76, 0.60), (0.79, 0.36))
    # ]:
    #     ax.annotate("", xy=b0, xytext=a0, arrowprops=arr)

    # # 虚线箭头改走框体上方，避免穿过框内
    # for x0 in (0.105, 0.295, 0.485):
    #     ax.annotate("", xy=(0.60, 0.37), xytext=(x0, 0.37),
    #                 arrowprops=dict(arrowstyle="->", lw=0.9, ls="--", color="0.45"))

    # ax.text(0.5, 0.01, "dashed arrows: associations tested in this study",
    #         ha="center", fontsize=7, color="0.35")

    # ax = fig.add_subplot(gs[1, 0])
    # y = np.arange(len(funnel))[::-1]
    # ax.barh(y, funnel["n"], color="0.55", height=0.62)
    # for yi, nv in zip(y, funnel["n"]):
    #     ax.text(nv + funnel["n"].max() * 0.02, yi, f"{int(nv)}", va="center", fontsize=7)
    # ax.set_yticks(y, [textwrap.fill(s, 20) for s in funnel["step"]])
    # ax.set_xlabel("Number of spectra"); ax.set_xlim(0, funnel["n"].max() * 1.2)
    # ax.set_title("QC funnel", loc="left"); panel(ax, "b")
    # ax = fig.add_subplot(gs[1, 1])
    # z = (D[tg] - D[tg].mean()) / D[tg].std(); z["family"] = D["family"]
    # prof = z.groupby("family").median(numeric_only=True); cnt = D["family"].value_counts(); prof = prof.loc[cnt.index]

    # # draw_heat(ax, prof.values, None, [tshort(t) for t in tg], [f"{fshort(f)} ({cnt[f]})" for f in prof.index],
    # #         vmin=-1.2, vmax=1.2, fmt="{:.1f}", fs=5.8, title="Median spectral profile (z-score)")
    
    # draw_heat(ax, prof.values, None, [tshort(t) for t in tg], [f"{fshort(f)} ({cnt[f]})" for f in prof.index],
    #         vmin=-1.2, vmax=1.2, fmt="{:.1f}", fs=5.8)
        
    # ax.set_title("Median spectral profile (z-score)", loc="left");
    # panel(ax, "c")
    # _save_pubfig(fig, out, "Figure1_framework_and_cohort")


    # ---------------- Figure 2: group-theory ceiling (no eta) ----------------
    fig, axs = plt.subplots(1, 3, figsize=(W, 2.9), layout="constrained")
    gtp = out / "tables" / "table_group_theory_ceiling.csv"
    ax = axs[0]
    if gtp.exists():
        g = pd.read_csv(gtp).dropna(subset=["n_raman_active_modes", "n_peak"]); g = g[g.n_raman_active_modes > 0]
        ax.scatter(g.n_raman_active_modes, g.n_peak, s=7, alpha=0.4, lw=0, color="0.35")
        lim = max(float(g.n_raman_active_modes.max()), float(g.n_peak.max()))
        ax.plot([0, lim], [0, lim], "k--", lw=0.8)
        if rho_gt is not None:
            annotate_stats(ax, f"Spearman $\\rho$ = {rho_gt:.2f}\nn = {len(g)}")
        ax.set_xlabel("Group-theoretical $N_{Raman}$"); ax.set_ylabel("Observed $N_{peak}$")
    ax.set_title("Mode availability vs peaks", loc="left"); panel(ax, "a")
    res = _peak_residual(D)
    D2 = D.copy()
    if res is not None:
        D2["peak_resid"] = res
    for ax, x, letter, ttl in ((axs[1], "raman_mode_fraction", "b", "Active-mode fraction"),
                               (axs[2], "degeneracy_fraction", "c", "Degeneracy fraction")):
        if res is None or x not in D2:
            ax.axis("off"); continue
        T = corr_table(D2, [x], ["peak_resid"], [c for c in COVARS if c in D2])
        _scatter_binned(ax, D2, T.rename(columns={}), x, "peak_resid")
        ax.set_ylabel("ln $N_{peak}$ residual | ln $N_{Raman}$")
        ax.set_title(ttl, loc="left"); panel(ax, letter)
    _save_pubfig(fig, out, "Figure2_group_theory_ceiling")

    # ---------------- Figures 3-4: bond geometry / chemistry ----------------
    for name, pairs, size in (("Figure3_bond_geometry", [("mean_reduced_mass", "n_peak"), ("mean_reduced_mass", "frac_low"),
                                                         ("mean_reduced_mass", "frac_high"), ("mean_delta_chi_bond", "w1_distance")], (W, 5.4)),
                              ("Figure4_chemical_heterogeneity", [("mixing_entropy", "gamma"), ("sigma_chi", "frac_high"),
                                                                  ("n_anion_group_types", "n_peak")], (W, 2.8))):
        nc = 2 if len(pairs) == 4 else 3
        fig, axs = plt.subplots(int(np.ceil(len(pairs) / nc)), nc, figsize=size, layout="constrained", squeeze=False)
        for ax, (x, y), letter in zip(axs.ravel(), pairs, "abcd"):
            _scatter_binned(ax, D, par, x, y); panel(ax, letter, dx=-0.16)
        for ax in axs.ravel()[len(pairs):]:
            ax.axis("off")
        _save_pubfig(fig, out, name)

    # ---------------- Figure 5: grouped-CV decomposition ----------------
    if Rt is not None and len(Rt):
        kt = [t for t in ["n_peak", "gamma", "domega_median", "r_overlap", "frac_low", "frac_mid", "frac_high", "w1_distance"]
              if t in set(Rt.target)]
        fig, ax = plt.subplots(figsize=(W, 3.6), layout="constrained")
        rows = [(full, short) for full, short in SET_SHORT if full in set(Rt.feature_set)]
        M = np.array([Rt[Rt.feature_set == full].set_index("target").reindex(kt).r2.to_numpy(float) for full, _ in rows])
        lim = max(0.3, float(np.nanmax(np.abs(M))))


        # im = draw_heat(ax, M, None, [tshort(t) for t in kt], [s for _, s in rows], vmin=-lim, vmax=lim, fmt="{:.2f}", fs=6.4,
        #                title="Grouped-CV $R^2$ by predictor set")

        im = draw_heat(ax, M, None, [tshort(t) for t in kt], [s for _, s in rows], vmin=-lim, vmax=lim, fmt="{:.2f}", fs=6.4,
                       )
        
        fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02).set_label("out-of-fold $R^2$")
        # panel(ax, "a")
        _save_pubfig(fig, out, "Figure5_grouped_cv_decomposition")

    # ---------------- Figure 6: joint group importance ----------------
    if Gdf is not None and len(Gdf):
        kt = [t for t in ["n_peak", "gamma", "domega_median", "r_overlap", "frac_low", "frac_mid", "frac_high", "w1_distance"] if t in Gdf.columns]
        rws = [GROUP_LABEL[g] for g in GROUPS if GROUP_LABEL[g] in Gdf.index]
        A = Gdf.reindex(rws)[kt]
        fig, ax = plt.subplots(figsize=(W, 3.3), layout="constrained")

        # im = draw_heat(ax, A.clip(lower=0).to_numpy(float), None, [tshort(t) for t in kt], rws,
        #                [GROUP_COLOR[g] for g in GROUPS if GROUP_LABEL[g] in rws], vmin=0,
        #                vmax=max(0.1, float(np.nanmax(A.to_numpy()))), cmap="YlOrRd", fmt="{:.2f}", annot_min=0.005, fs=6.4,
        #                title="Joint held-out permutation importance")
        

        im = draw_heat(ax, A.clip(lower=0).to_numpy(float), None, [tshort(t) for t in kt], rws,
                       [GROUP_COLOR[g] for g in GROUPS if GROUP_LABEL[g] in rws], vmin=0,
                       vmax=max(0.1, float(np.nanmax(A.to_numpy()))), cmap="YlOrRd", fmt="{:.2f}", annot_min=0.005, fs=6.4,
                       )

        fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02).set_label("$\\Delta R^2$")
        # panel(ax, "a")
        _save_pubfig(fig, out, "Figure6_group_importance")

    # ---------------- Figure 7: conceptual framework ----------------
    fig = plt.figure(figsize=(W, 4.3)); ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
    stages = [(0.04, 0.56, 0.26, 0.24, "Crystal structure", "size / complexity\nsymmetry / group theory\nbond geometry"),
              (0.37, 0.56, 0.26, 0.24, "Composition-derived chemistry", "elemental composition\nsite mixing\nelectronegativity"),
              (0.70, 0.56, 0.26, 0.24, "Raman mode space", "allowed modes\nfrequency distribution\nmode degeneracy"),
              (0.22, 0.14, 0.26, 0.24, "Spectral complexity", "peak count\nspacing\nband fractions"),
              (0.57, 0.14, 0.26, 0.24, "Band resolvability", "linewidth\noverlap\nresolved peaks")]
    for x, y, w, h, title, body in stages:
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=.012", fill=False, lw=1.0, ec="0.25"))
        ax.text(x + w / 2, y + h - 0.055, title, ha="center", va="center", fontsize=9, fontweight="bold")
        ax.text(x + w / 2, y + h / 2 - 0.01, body, ha="center", va="center", fontsize=7.4, linespacing=1.45)
    for a, b in [((.31, .68), (.36, .68)), ((.64, .68), (.69, .68)), ((.83, .55), (.69, .39)), ((.45, .55), (.38, .39)), ((.49, .265), (.56, .265))]:
        ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="->", lw=1.15))
    # ax.text(0.5, 0.04, "Proposed framework (hypothesis): arrows are tested only via the associations in Figures 2-6",
            # ha="center", va="center", fontsize=8.2, fontweight="bold")
    _save_pubfig(fig, out, "Figure7_integrated_framework", layout=False)

    # ---------------- Supplementary ----------------
    fig = plt.figure(figsize=(W, 3.4), layout="constrained")
    gs = fig.add_gridspec(1, 2, width_ratios=[0.85, 1.25])
    ax = fig.add_subplot(gs[0]); y = np.arange(len(funnel))[::-1]
    ax.barh(y, funnel.n, color="0.55", height=0.6)
    for yi, nv in zip(y, funnel.n):
        ax.text(nv + funnel.n.max() * 0.02, yi, str(int(nv)), va="center", fontsize=7)
    ax.set_yticks(y, [textwrap.fill(s, 20) for s in funnel.step]); ax.set_xlabel("Number of spectra"); ax.set_xlim(0, funnel.n.max() * 1.2)
    ax.set_title("QC funnel", loc="left"); panel(ax, "a")
    ax = fig.add_subplot(gs[1])
    draw_heat(ax, prof.values, None, [tshort(t) for t in tg], [f"{fshort(f)} ({cnt[f]})" for f in prof.index],
              vmin=-1.2, vmax=1.2, fmt="{:.1f}", fs=5.8, title="Median spectral profile (z-score)")
    panel(ax, "b")
    _save_suppfig(fig, out, "FigureS1_qc_and_family_profile")

    names = [f.name for f in FEATURES if f.name in D and D[f.name].notna().sum() >= 30 and D[f.name].nunique() > 1]
    order, cols, seps = feat_rows(names)                 # 传入列名(不是 Feat 对象) -> 修复 S2 全空
    fig, ax = plt.subplots(figsize=(W, 8.4), layout="constrained")
    im = draw_heat(ax, to_mat(par, order, tg, "rho"), to_mat(par, order, tg, "q"), [tshort(t) for t in tg],
                   [flabel(c) for c in order], cols, vmin=-0.6, vmax=0.6, fmt="{:.2f}", sep_rows=seps, fs=5.5,
                   title="Partial Spearman correlations (adj. ln N_atom, ln SNR)")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02).set_label("partial $\\rho$")
    ax.set_xlabel("* q<0.05, ** q<0.01, *** q<0.001 (BH-FDR over all cells)", fontsize=7)
    panel(ax, "a")
    _save_suppfig(fig, out, "FigureS2_full_partial_correlations")

    if cat is not None and len(cat):
        facs = list(cat.factor.dropna().unique())
        M = cat.pivot(index="factor", columns="target", values="eta2_H").reindex(index=facs, columns=tg)
        Q = cat.pivot(index="factor", columns="target", values="q").reindex(index=facs, columns=tg)
        fig, ax = plt.subplots(figsize=(W, 3.2), layout="constrained")
        im = draw_heat(ax, M.to_numpy(float), Q.to_numpy(float), [tshort(t) for t in tg], list(M.index), vmin=0,
                       vmax=max(0.3, float(np.nanmax(M.to_numpy()))), cmap="Purples", fmt="{:.2f}", fs=6.2,
                       title="Kruskal-Wallis $\\eta^2_H$")
        fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02).set_label("$\\eta^2_H$")
        panel(ax, "a")
        _save_suppfig(fig, out, "FigureS3_categorical_effects")

    if S is not None and len(S):
        resp = [y for y in tg if y in set(S.response)]
        if resp:
            fig = _forest_grid(S, resp)
            _save_suppfig(fig, out, "FigureS4_regression_forest")

    if L is not None and len(L):
        A = L[L.feature_set == "all"]
        piv = A.pivot(index="held_out_family", columns="target", values="skill_vs_train_mean")
        fig, ax = plt.subplots(figsize=(W, 4.2), layout="constrained")
        im = draw_heat(ax, piv.to_numpy(float), None, [tshort(t) for t in piv.columns], [fshort(f) for f in piv.index],
                       vmin=-1, vmax=1, fmt="{:.2f}", fs=6.0, title="Leave-one-family-out skill score")
        fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02).set_label("skill score")
        panel(ax, "a")
        _save_suppfig(fig, out, "FigureS5_lofo_extrapolation")

    fig = plt.figure(figsize=(W, 7.4), layout="constrained")
    gs = fig.add_gridspec(3, 1, height_ratios=[1.7, 0.9, 0.9])
    ax = fig.add_subplot(gs[0])
    if C is not None and not C.empty:
        show = C.loc[C.abs().max(axis=1).sort_values(ascending=False).head(18).index]
        lim = max(0.2, float(show.abs().max().max()))
        im = draw_heat(ax, show.to_numpy(float), None, [tshort(t) for t in show.columns],
                       [("Missingness: " + flabel(x[3:])) if x.startswith("NA:") else flabel(x) for x in show.index],
                       vmin=-lim, vmax=lim, fmt="{:+.2f}", fs=5.7, title="Lasso coefficients (top 18)")
        fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02).set_label("standardised coefficient")
    panel(ax, "a")
    ax = fig.add_subplot(gs[1])
    if q2:
        ax.plot(np.arange(1, len(q2) + 1), q2, "o-", lw=1.2, ms=3.5)
        ax.set_xlabel("PLS latent variables"); ax.set_ylabel("Cumulative $Q^2$")
        ax.set_xticks(range(1, len(q2) + 1)); ax.axhline(0, color="0.6", lw=0.6)
    ax.set_title("PLS grouped-CV $Q^2$", loc="left"); panel(ax, "b")
    ax = fig.add_subplot(gs[2])
    if sens is not None and len(sens):
        ss = sens.set_index("spec"); vals = ss["spearman_of_all_rho_vs_main"].to_numpy(float)
        ax.bar(np.arange(len(vals)), vals, color="0.55"); ax.axhline(1, color="0.35", ls="--", lw=0.7)
        ax.set_xticks(np.arange(len(vals)), list(ss.index)); ax.set_ylim(0, 1.05)
        ax.set_ylabel("Rank agreement with main")
    else:
        ax.text(0.5, 0.5, "QC sensitivity not run (--fast / --skip-sensitivity)", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
    ax.set_title("QC-threshold sensitivity", loc="left"); panel(ax, "c")
    _save_suppfig(fig, out, "FigureS6_sparse_pls_sensitivity")

    cap = ["# Figure captions (descriptive; result statements are in paper_results.md)\n"]
    cap += [f"## {k}\n{v}\n" for k, v in {**MAIN_CAPTIONS, **SUPP_CAPTIONS}.items()]
    (out / "figure_captions.md").write_text("\n".join(cap), encoding="utf-8")
    man = {"main": sorted(p.stem for p in (out / "figures" / "main").glob("*.png")),
           "supplementary": sorted(p.stem for p in (out / "figures" / "supplementary").glob("*.png"))}
    (out / "figure_manifest.json").write_text(json.dumps(man, indent=2), encoding="utf-8")


def write_summary(out, args, funnel, n, gcol, par, cat, S, Rt, Dt, G, lasso, Fdf, Gdf, q2, sens, rho_gt, tg):
    L = [f"# Auto-generated results summary\n\nAll numbers below come from this run (seed={args.seed}); interpret with domain knowledge.\n",
         f"## 0. Cohort\nn = {n} spectra after QC; grouping unit for CV / clustering: `{gcol}`.\n", md_table(funnel), ""]

    # Data-driven statements only (no pre-written mechanism).
    L += ["## Key findings (auto-generated from this run's numbers)\n"]
    hp = _headline_pairs(par)
    if len(hp):
        L += ["Hypothesis-driven pairs fixed in code (HEADLINE_PAIRS); report as exploratory if chosen after viewing results.\n",
              md_table(pd.DataFrame({"structure": hp.feature.map(flabel), "spectrum": hp.target.map(tlabel),
                                     "partial rho [95% CI]": [_fmt_ci(r) for r in hp.itertuples()],
                                     "q": hp.q.map(fmt_p)})), ""]
    L += [f"- {s}" for s in _evidence_sentences(par, Rt, Dt, Gdf, rho_gt, sens)] + [""]
    top = par.dropna(subset=["rho"]).assign(a=lambda d: d.rho.abs()).sort_values("a", ascending=False).head(15)
    L += ["## 1. Strongest partial correlations (adjusted for ln N_atom, ln SNR)\n",
          md_table(pd.DataFrame({"structure": top.feature.map(flabel), "spectrum": top.target.map(tlabel),
                                 "partial rho [95% CI]": [f"{r.rho:+.2f} [{r.ci_lo:+.2f}, {r.ci_hi:+.2f}]" for r in top.itertuples()],
                                 "q": top.q.map(fmt_p), "n": top.n.astype(int)})), ""]
    g = cat.dropna(subset=["eta2_H"]).sort_values("eta2_H", ascending=False).head(10)
    L += ["## 2. Categorical factors (top 10 eta2_H)\n",
          md_table(pd.DataFrame({"factor": g.factor, "spectrum": g.target.map(tlabel), "eta2_H": g.eta2_H.round(2), "q": g.q.map(fmt_p)})), ""]
    if S is not None and len(S):
        s = S[(S.q < FDR_ALPHA) & S.tier.isin(["M2 extended", "M3 extended + family FE"])].copy()
        s["abs"] = s.beta_std.abs()
        s = s.sort_values("abs", ascending=False).head(20)
        L += ["## 3. Regression: significant standardised coefficients (BH q<0.05; M2 vs M3 tells whether an effect survives within-family)\n",
              md_table(pd.DataFrame({"tier": s.tier, "response": s.response.map(tlabel), "term": s.term.map(flabel),
                                     "beta [95% CI]": [f"{r.beta_std:+.2f} [{r.ci_lo:+.2f}, {r.ci_hi:+.2f}]" for r in s.itertuples()],
                                     "q": s.q.map(fmt_p), "R2": s.r2.round(2)})), ""]
    P = Rt[~Rt.feature_set.str.startswith("All minus")].copy()
    P["cell"] = [f"{r.r2:.2f} [{r.ci_lo:.2f}, {r.ci_hi:.2f}]" for r in P.itertuples()]
    Pv = P.pivot(index="target", columns="feature_set", values="cell").reindex(index=tg, columns=list(dict.fromkeys(P.feature_set)))
    Pv.index = [tlabel(t) for t in Pv.index]
    L += ["## 4. Out-of-fold R2 [95% cluster-bootstrap CI], grouped CV\n", md_table(Pv.reset_index().rename(columns={"target": "spectrum", "index": "spectrum"})), ""]
    if Dt is not None and len(Dt):
        d = Dt[Dt.contrast.isin(["All - Categorical", "All - Measurement"])].copy()
        d["cell"] = [f"{r.delta_r2:+.2f} [{r.ci_lo:+.2f}, {r.ci_hi:+.2f}]" for r in d.itertuples()]
        L += ["### Added value of structure/chemistry over baselines (paired dR2)\n",
              md_table(d.pivot(index="target", columns="contrast", values="cell").reindex(tg).rename(index=tlabel).rename_axis("spectrum").reset_index()), ""]
    if G is not None and len(G):
        gg = G[G.feature_set == "all"].copy()
        gg["generalisation gap"] = gg["grouped_cv_r2_same_rows"] - gg["lofo_r2_pooled"]
        L += ["### Leave-one-family-out extrapolation (all features)\n",
              md_table(gg[["target", "lofo_skill_pooled", "lofo_r2_pooled", "grouped_cv_r2_same_rows", "generalisation gap"]]
                       .assign(target=lambda x: x.target.map(tlabel))), ""]
    if lasso:
        L += ["## 5. Sparse relations (nested-CV Lasso; percentages = bootstrap selection frequency)\n"]
        L += [f"- **{tlabel(t)}** (R2_cv={o['r2']:.2f} [{o['ci'][0]:.2f}, {o['ci'][1]:.2f}], {o['nnz']} terms): `{o['eq']}`" for t, o in lasso.items()]
        L.append("")
    if Gdf is not None:
        Gd = Gdf.copy()
        Gd.columns = [tlabel(t) for t in Gd.columns]
        L += ["## 6. Joint group-permutation importance (dR2)\n", md_table(Gd.round(3).reset_index().rename(columns={"index": "group"})), ""]
    L += ["## 7. PLS cumulative Q2 (grouped CV)\n" + ", ".join(f"{k + 1} LV = {v:.3f}" for k, v in enumerate(q2)), ""]
    if rho_gt is not None:
        L += [f"## 8. Group-theory ceiling\nSpearman rho(N_Raman, N_peak) = {rho_gt:.2f}.\n"]
    if sens is not None and len(sens):
        L += ["## 9. Sensitivity of partial correlations to QC thresholds\n", md_table(sens), ""]
    (out / "results_summary.md").write_text("\n".join(L), encoding="utf-8")


def write_methods(out, args, n, gcol, funnel, n_feats):
    txt = f"""# Statistical methods — publication-oriented workflow

**Study design and QC.** The analysis began with the complete RRUFF-derived cohort and retained spectra with quality in {{{', '.join(args.qualities)}}}, SNR >= {args.min_snr:g}, at least {args.min_peaks} detected peaks, and {args.orientation} orientation (final n = {n}). The grouping unit for all resampling-based inference and prediction was chemical formula, preventing replicate spectra of the same composition from being split across folds. Because the minimum-peak criterion conditions on a response-related quantity, sensitivity analyses repeated the association screen at SNR thresholds of 10 and 40 and peak-count thresholds of 2 and 5.

**Mechanistic feature hierarchy.** Structural and composition-derived descriptors were explicitly partitioned into size/complexity, symmetry/group theory, bond geometry, composition-derived chemistry, and local environment (SOAP). Group-theoretical descriptors were retained only when the reported crystallographic symmetry was internally consistent; otherwise the corresponding symmetry block was masked. Spectral endpoints represented complementary aspects of band complexity and resolvability: peak count, intensity evenness, mean linewidth, median peak spacing, peak overlap, low/mid/high spectral fractions, and Wasserstein distance. Mapping efficiency eta = N_peak/N_Raman was used only as an inferential endpoint and excluded from machine-learning predictors to avoid definitional leakage.

**Primary association analysis.** Rank-based partial correlations adjusted for ln N_atom and ln SNR were used to separate structural effects from size and measurement quality. Confidence intervals were calculated on the rank scale and p-values were controlled by Benjamini–Hochberg FDR at {FDR_ALPHA}. Categorical effects were quantified with Kruskal–Wallis eta-squared. The primary interpretation emphasizes effect size and confidence intervals rather than significance alone.

**Multivariable inference.** Standardised OLS models were fitted as a core model (M1), an extended mechanistic model (M2), and an extended model with chemical-family fixed effects (M3). Chemical-formula cluster-robust standard errors were used for inference. Persistence of an effect from M2 to M3 was used as a within-family robustness check, not as evidence of causality. VIF diagnostics were retained to monitor collinearity.

**Predictive validation and incremental value.** HistGradientBoosting regressors were evaluated using repeated {args.cv_folds}-fold grouped CV. All preprocessing was learned within training folds only. Out-of-fold R2 and paired R2 differences were quantified with chemical-formula cluster bootstrap confidence intervals. Two baselines were prespecified: measurement conditions (SNR and excitation wavelength) and categorical descriptors (chemical family and crystal system). The key publication comparison is the incremental dR2 of all structure/chemistry descriptors over the categorical baseline. Leave-one-family-out validation was used as an extrapolation stress test and is reported as a secondary generalisation analysis.

**Feature-group attribution.** Single-feature and joint group permutation importance were computed on held-out folds. Joint permutation preserves within-group feature covariance and asks whether an entire mechanistic descriptor class contributes out-of-sample information. Drop-group ablation provides a complementary dR2 estimate.

**Sparse mechanistic models.** Nested grouped-CV Lasso used the one-standard-error rule to obtain compact relations. Selection frequency was estimated by chemical-formula cluster bootstrap. Missingness indicators are labelled explicitly as “Missingness: ...” and are interpreted as data-availability signals rather than physical descriptors.

**Latent structure and robustness.** PLS cumulative Q2 was estimated under grouped CV, with imputation/standardisation fitted within each training fold. QC sensitivity was used as a robustness check. UMAP, when requested, is exploratory only and is not used for inferential conclusions.
"""
    (out / "methods_text.md").write_text(txt, encoding="utf-8")


# =============================================================================
# 11. main
# =============================================================================
def main(argv=None) -> int:
    global DPI, FIG_TITLES
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", type=Path, default=Path("results/dataset.csv"))
    ap.add_argument("--out", type=Path, default=Path("results"))
    ap.add_argument("--qualities", nargs="+", default=["excellent", "fair"])
    ap.add_argument("--min-snr", type=float, default=20)
    ap.add_argument("--min-peaks", type=int, default=3)
    ap.add_argument("--orientation", default="unoriented", help="'any' 表示不按取向过滤")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--cv-folds", type=int, default=5)
    ap.add_argument("--repeats", type=int, default=3, help="重复分组 CV 次数")
    ap.add_argument("--n-boot", type=int, default=300, help="聚类 bootstrap 次数 (CI)")
    ap.add_argument("--n-perm", type=int, default=5, help="置换重要性重复次数")
    ap.add_argument("--n-stab", type=int, default=100, help="Lasso 稳定性选择 bootstrap 次数")
    ap.add_argument("--min-family-size", type=int, default=15)
    ap.add_argument("--fast", action="store_true", help="调试: 1 次重复, 少量 bootstrap, 跳过消融/LOFO/稳健性")
    ap.add_argument("--skip-sensitivity", action="store_true")
    ap.add_argument("--dpi", type=int, default=300, help="PNG 分辨率")
    ap.add_argument("--fig-titles", action="store_true", help="在图内显示总标题 (默认不显示, 交给图注)")
    ap.add_argument("--paper-mode", action="store_true", default=True, help="生成论文核心综合图/表/结果摘要（默认开启）")
    args = ap.parse_args(argv)
    DPI, FIG_TITLES = args.dpi, args.fig_titles
    if args.fast:
        args.repeats, args.n_boot, args.n_perm, args.n_stab = 1, 100, 2, 30

    out = args.out
    tab, rep = out / "tables", out / "regression_reports"
    for p in (tab, rep):
        p.mkdir(parents=True, exist_ok=True)
    set_style()
    # Analysis steps may construct transient diagnostic figures for in-run checks, but never save them.
    # Only publication_figure_pipeline() writes PNG files: Figure 1–7 and Supplementary Figure S1–S6.
    save = lambda fig, name: plt.close(fig)

    raw_df = read_raw(args.csv)
    d, funnel = apply_qc(raw_df, args.qualities, args.min_snr, args.min_peaks, args.orientation)
    print(funnel.to_string(index=False))
    if len(d) < 100:
        print(f"[warning] only {len(d)} rows survive QC; many analyses will be unstable")
    R, D, rec = build_frames(d)
    gcol = group_column(D)
    feats = [f.name for f in FEATURES if f.name in D and D[f.name].notna().sum() >= 30 and D[f.name].nunique() > 1]
    tg = [t.name for t in TARGETS if t.name in D and D[t.name].notna().sum() >= 30]
    absent = [f.name for f in FEATURES if f.name not in feats]
    if absent:
        print(f"[warning] descriptors absent/constant and skipped: {absent}")
    covars = [c for c in COVARS if c in D]
    rec.to_csv(tab / "table_transforms.csv", index=False)
    D.to_csv(tab / "analysis_frame_transformed.csv", index=False)
    print(f"n = {len(D)}, {len(feats)} descriptors, {len(tg)} spectral targets, grouping = {gcol}")
    print(f"  symmetry_consistent=False: {int((~R['symmetry_consistent']).sum())} rows (group-theory columns masked)")

    print("[S0] cohort");                 step0_cohort(R, D, funnel, tg, save, tab)
    print("[S1] correlations");           raw, par = step1_correlation(D, feats, tg, save, tab, covars)
    print("[S2] categorical effects");    cat = step2_categorical(D, tg, save, tab)
    print("[S3] key relationships");      step3_key_scatter(D, par, raw, feats, save, tab)
    print("[S4] group-theory ceiling");   rho_gt = step4_group_theory(R, save, tab)
    print("[S5] regression");             S = step5_regression(D, tg, gcol, save, tab, rep)
    mtg = [t for t in tg if t not in ML_EXCLUDE]            # 机器学习环节的目标 (不含 eta)
    print("[S6] predictive power (grouped CV)")
    Rt, Dt, oof_all = step6_predictive(D, feats, mtg, gcol, args, save, tab)
    print("[S6b] leave-one-family-out");  L, G = step6b_lofo(D, feats, mtg, gcol, args, oof_all, tab)
    print("[S7] permutation importance"); Fdf, Gdf = step7_importance(D, feats, mtg, gcol, args, save, tab, Dt)
    print("[S8] Lasso (nested CV + stability)"); lasso, C, Fq = step8_lasso(D, feats, mtg, gcol, args, save, tab)
    print("[S9] PLS");                    q2 = step9_pls(D, feats, mtg, gcol, args, save, tab)
    sens = None
    if not (args.fast or args.skip_sensitivity):
        print("[S10] QC sensitivity");    sens = step10_sensitivity(raw_df, args, feats, tg, covars, tab)

    # 文本产出先写 (不依赖绘图), 再出投稿图: 出图失败不会丢失 summary / methods
    write_summary(out, args, funnel, len(D), gcol, par, cat, S, Rt, Dt, G, lasso, Fdf, Gdf, q2, sens, rho_gt, mtg)
    write_methods(out, args, len(D), gcol, funnel, len(feats))
    write_paper_results(out, args, len(D), par, S, Rt, Dt, Gdf, rho_gt, sens)
    print("[PAPER] publication figure pipeline: 7 main + 6 supplementary (PNG only)")


    try:
        publication_figure_pipeline(out, args, R, D, funnel, par, cat, S, Rt, Dt, L, G, lasso, C, Fq, q2, Gdf, rho_gt, sens, mtg)
    except Exception:
        import traceback
        traceback.print_exc()
        print("[error] publication figure pipeline failed; analysis tables / summary were already written")
    qa_lines = [f"{n}: {k} overlapping text pairs" + ("".join(f"\n    - {i}" for i in iss) if iss else "") for n, k, iss in QA_LOG]
    (out / "figure_qa.txt").write_text("\n".join(qa_lines) + "\n", encoding="utf-8")
    n_bad = sum(1 for _, k, iss in QA_LOG if iss)
    print(f"[figure QA] {len(QA_LOG)} figures checked, {n_bad} with issues -> {out / 'figure_qa.txt'}")
    (out / "run_config.json").write_text(json.dumps(
        dict(args={k: str(v) for k, v in vars(args).items()}, python=sys.version.split()[0], numpy=np.__version__,
             pandas=pd.__version__, n=len(D)), indent=2), encoding="utf-8")
    print(f"\n[paper] main figures      : {out / 'figures' / 'main'}")
    print(f"[paper] supplementary figs : {out / 'figures' / 'supplementary'}")
    print(f"[paper] captions           : {out / 'figure_captions.md'}")
    print(f"[paper] key table   : {out / 'tables' / 'table_paper_key_results.csv'}")
    print(f"[paper] narrative   : {out / 'paper_results.md'}")
    print(f"\n[done] outputs in {out.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
