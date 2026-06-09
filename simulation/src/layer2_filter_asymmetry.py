"""
===============================================================================
Layer 2 -- Maternal filtering and directional reproductive isolation
===============================================================================
Project : ovule-gated pollen competition and passenger RI
Builds on: layer1_competition_engine.py (the O -> K hinge)
Maps to : theory note sections "Maternal filtering", "Prezygotic reproductive
          isolation", "Directional asymmetry", and simulation Layers 2 + 5.

WHAT THIS MODULE ADDS TO LAYER 1
--------------------------------
Layer 1 established the competition ratio K = P/O and the selection it creates.
Layer 2 puts a *maternal filter* on top and turns it into directional RI:

  1. evolved filter strictness        F*(K) = Fmax * K / (K + K50)     (rises with K)
  2. pollen-pistil mismatch           D_ij  = delta_D * G_ij           (grows with divergence)
  3. directional RI, PRIMARY index    R_{i<-j} = tanh(F_i * D_ij / 2)  (Sobel-Chen RI4)
     directional RI, SENSITIVITY      R_alt    = 1 - exp(-F_i * D_ij)  (1 - H/C)
  4. reciprocal asymmetry             A_ij = R_{i<-j} - R_{j<-i}

The maternal filter is the only *directional* ingredient: D_ij is shared by a
pair, so the two cross directions can differ only through F_i vs F_j.

PRIMARY vs SENSITIVITY INDEX
----------------------------
Per the project decision, the Sobel-Chen directional index R = tanh(F D / 2) is
PRIMARY. The simpler R = 1 - H/C is retained as a SENSITIVITY check. They share
the sign theorem and the hump-shape but differ in scale: at small mismatch the
Sobel-Chen asymmetry carries a factor 1/2 (A ~ 1/2 * D * (F_i - F_j)) while the
1 - H/C form linearises to A ~ D * (F_i - F_j) -- twice the slope.

WHAT IS VALIDATED HERE (assertions halt the run on failure)
-----------------------------------------------------------
  (V1) sign theorem      : sign(A_ij) = sign(F_i - F_j) exactly, both indices.
  (V2) linearisation     : A_SC/[D dF] -> 1/2 ; A_alt/[D dF] -> 1 as D -> 0.
  (V3) hump-shape        : A(D) -> 0 at D -> 0 and D -> infinity; one interior peak.
  (V4) interaction form  : simulating pairs and regressing A on its candidate
                           predictors recovers G_ij*(logO_j - logO_i) as the
                           dominant positive term -- the signature the real-data
                           discriminating regression will look for.

Genetic architectures (independent / linked / pleiotropic) are Layer 3 and are
deliberately NOT here; this module fixes D_ij = delta_D * G_ij to leading order.
===============================================================================
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Reuse the validated Layer 1 competition ratio; fall back to a local copy so
# this file also runs standalone.
try:
    from layer1_competition_engine import K_from_O
except Exception:  # pragma: no cover
    def K_from_O(O, psi, c, alpha):
        O = np.asarray(O, dtype=float)
        return psi * c * O ** (alpha - 1.0)


# =============================================================================
# Core deterministic functions
# =============================================================================
def Fstar(K: np.ndarray, Fmax: float = 2.0, K50: float = 2.0) -> np.ndarray:
    """Evolved maternal filter strictness, F*(K) = Fmax K / (K + K50).

    A Hill/saturating function: zero when K = 0 (no surplus pollen, no point
    being choosy), rising monotonically toward Fmax as K grows. dF*/dK > 0 for
    all K >= 0, which is the only property the directional prediction needs.
    """
    K = np.asarray(K, dtype=float)
    return Fmax * K / (K + K50)


def R_sobelchen(F: np.ndarray, D: np.ndarray) -> np.ndarray:
    """PRIMARY directional RI index: Sobel-Chen RI4 = tanh(F D / 2).

    Identical to (1 - h)/(1 + h) with h = exp(-F D); bounded in [0, 1) for
    F, D >= 0; equals 0 when there is no mismatch and approaches 1 as the
    maternal filter x mismatch product grows.
    """
    return np.tanh(np.asarray(F) * np.asarray(D) / 2.0)


def R_oneminus_hc(F: np.ndarray, D: np.ndarray) -> np.ndarray:
    """SENSITIVITY directional RI index: 1 - H/C = 1 - exp(-F D)."""
    return 1.0 - np.exp(-np.asarray(F) * np.asarray(D))


def mismatch(G: np.ndarray, delta_D: float) -> np.ndarray:
    """Leading-order pollen-pistil mismatch, D_ij = delta_D * G_ij (symmetric)."""
    return delta_D * np.asarray(G, dtype=float)


def asymmetry(Rfun, Fi: np.ndarray, Fj: np.ndarray, D: np.ndarray) -> np.ndarray:
    """Reciprocal asymmetry A_ij = R(F_i, D) - R(F_j, D) for a given index."""
    return Rfun(Fi, D) - Rfun(Fj, D)


# =============================================================================
# Validation block (the heart of Layer 2)
# =============================================================================
def validate_core(seed: int = 0) -> None:
    rng = np.random.default_rng(seed)

    # V1 -- sign theorem, both indices, over random (F_i, F_j, D>0)
    Fi = rng.uniform(0, 4, 200_000); Fj = rng.uniform(0, 4, 200_000)
    D = rng.uniform(1e-3, 4, 200_000)
    for name, Rf in [("SobelChen", R_sobelchen), ("1-H/C", R_oneminus_hc)]:
        A = asymmetry(Rf, Fi, Fj, D)
        # sign(A) must equal sign(Fi - Fj); ties (Fi==Fj) give A==0, also fine
        ok = np.all(np.sign(A) == np.sign(Fi - Fj))
        assert ok, f"V1 sign theorem FAILED for {name}"

    # V2 -- small-D linearisation slopes: 1/2 (Sobel-Chen), 1 (1 - H/C)
    Fi0, Fj0, Dsm = 1.5, 0.9, 1e-3
    s_sc = asymmetry(R_sobelchen, Fi0, Fj0, Dsm) / (Dsm * (Fi0 - Fj0))
    s_alt = asymmetry(R_oneminus_hc, Fi0, Fj0, Dsm) / (Dsm * (Fi0 - Fj0))
    assert abs(s_sc - 0.5) < 1e-3, f"V2 Sobel-Chen slope {s_sc} != 1/2"
    assert abs(s_alt - 1.0) < 1e-2, f"V2 1-H/C slope {s_alt} != 1"

    # V3 -- hump-shape: A(D) vanishes at both ends, single interior maximum
    Dg = np.linspace(1e-3, 50, 20_000)
    for name, Rf in [("SobelChen", R_sobelchen), ("1-H/C", R_oneminus_hc)]:
        A = asymmetry(Rf, 1.5, 0.9, Dg)
        k = int(A.argmax())
        assert 0 < k < len(Dg) - 1, f"V3 peak not interior for {name}"
        assert A[0] < A[k] and A[-1] < A[k], f"V3 not hump-shaped for {name}"
        assert A[-1] < 1e-3, f"V3 high-D asymmetry did not collapse for {name}"
    print("[validate] V1 sign theorem, V2 linearisation (1/2 vs 1), "
          "V3 hump-shape: PASSED for both indices.")


# =============================================================================
# Simulating lineage pairs and recovering the interaction (V4)
# =============================================================================
@dataclass
class Layer2Config:
    n_lineages: int = 120
    psi: float = 0.02
    c: float = 500.0
    alpha: float = 0.5          # < 1 so K falls with ovule number (Layer 1 hinge)
    Fmax: float = 2.0
    K50: float = 2.0
    delta_D: float = 5.0        # mismatch per unit genetic distance
    O_min: float = 2.0
    O_max: float = 20000.0
    G_mean: float = 0.15        # mean K2P-like divergence for simulated pairs
    n_pairs: int = 3000
    seed: int = 11


def simulate_pairs(cfg: Layer2Config) -> pd.DataFrame:
    """Draw lineages (ovule numbers -> K -> F) and random pairs, compute A_ij.

    Genetic distance G is drawn independently of the ovule contrast so the
    regression in V4 can separate the interaction G*dlogO from its main effects.
    """
    rng = np.random.default_rng(cfg.seed)
    # lineages: log-uniform ovule numbers spanning few..many
    O = np.exp(rng.uniform(np.log(cfg.O_min), np.log(cfg.O_max), cfg.n_lineages))
    K = K_from_O(O, cfg.psi, cfg.c, cfg.alpha)
    F = Fstar(K, cfg.Fmax, cfg.K50)

    i = rng.integers(0, cfg.n_lineages, cfg.n_pairs)
    j = rng.integers(0, cfg.n_lineages, cfg.n_pairs)
    keep = i != j
    i, j = i[keep], j[keep]
    G = np.abs(rng.normal(cfg.G_mean, cfg.G_mean * 0.5, i.size))  # >=0 divergence
    D = mismatch(G, cfg.delta_D)

    A_sc = asymmetry(R_sobelchen, F[i], F[j], D)
    A_alt = asymmetry(R_oneminus_hc, F[i], F[j], D)
    dlogO = np.log(O[j]) - np.log(O[i])          # positive when i is the few-ovule mother
    return pd.DataFrame(dict(
        i=i, j=j, O_i=O[i], O_j=O[j], K_i=K[i], K_j=K[j], F_i=F[i], F_j=F[j],
        G=G, D=D, dlogO=dlogO, GxdlogO=G * dlogO, A_sobelchen=A_sc, A_oneminus_hc=A_alt,
    ))


def _ols(y: np.ndarray, X: np.ndarray):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = len(y) - X.shape[1]
    s2 = resid @ resid / dof
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    return beta, beta / se


def validate_interaction(df: pd.DataFrame) -> pd.DataFrame:
    """V4: regress simulated A on [GxdlogO, dlogO, G]; the interaction must lead.

    The derived form is A_ij ~ beta1 * G_ij * (logO_j - logO_i). So with
    standardised predictors the GxdlogO coefficient should be positive and the
    largest in magnitude, while the main effects sit near zero.
    """
    y = df["A_sobelchen"].values
    cols = ["GxdlogO", "dlogO", "G"]
    Z = df[cols].values
    Zs = (Z - Z.mean(0)) / Z.std(0)                 # standardise for comparability
    ys = (y - y.mean()) / y.std()
    X = np.column_stack([np.ones(len(ys)), Zs])
    beta, t = _ols(ys, X)
    out = pd.DataFrame({"predictor": ["intercept"] + cols,
                        "std_coef": beta, "t": t})
    inter = out.loc[out.predictor == "GxdlogO"].iloc[0]
    biggest = out.iloc[1:].reindex(out.iloc[1:]["std_coef"].abs().sort_values().index)
    assert inter["std_coef"] > 0, "V4 interaction coefficient is not positive"
    assert biggest.iloc[-1]["predictor"] == "GxdlogO", \
        "V4 interaction is not the dominant predictor"
    print("[validate] V4 interaction recovery: G*dlogO is the dominant positive "
          f"predictor (std coef {inter['std_coef']:.3f}, t {inter['t']:.1f}).")
    return out


# =============================================================================
# Figures
# =============================================================================
def fig_hump(path: Path, Fi: float = 1.5, Fj: float = 0.9) -> None:
    D = np.linspace(0, 12, 600)
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.plot(D, asymmetry(R_sobelchen, Fi, Fj, D), color="black", lw=2,
            label="Sobel–Chen $\\tanh(FD/2)$  (primary)")
    ax.plot(D, asymmetry(R_oneminus_hc, Fi, Fj, D), color="#C0392B", lw=2, ls="--",
            label="$1-H/C$  (sensitivity)")
    ax.axhline(0, color="grey", lw=0.8)
    ax.set(xlabel="mismatch $D_{ij}$", ylabel="reciprocal asymmetry $A_{ij}$",
           title=f"Asymmetry is hump-shaped in divergence  ($F_i={Fi}>F_j={Fj}$)")
    ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def fig_interaction(df: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 5))
    qs = df["G"].quantile([0.15, 0.5, 0.85]).values
    cmap = plt.get_cmap("viridis")
    for k, g in enumerate(qs):
        sub = df[np.abs(df["G"] - g) < 0.02].sort_values("dlogO")
        if len(sub) > 5:
            ax.scatter(sub["dlogO"], sub["A_sobelchen"], s=8,
                       color=cmap(k / 2), label=f"G ≈ {g:.2f}")
    ax.axhline(0, color="grey", lw=0.8); ax.axvline(0, color="grey", lw=0.8)
    ax.set(xlabel="$\\log O_j - \\log O_i$  (positive: $i$ is the few-ovule mother)",
           ylabel="$A_{ij}$ (Sobel–Chen)",
           title="Asymmetry follows the ovule contrast, scaled by divergence")
    ax.legend(frameon=False, title="divergence slice")
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def fig_filter(cfg: Layer2Config, path: Path) -> None:
    O = np.logspace(np.log10(cfg.O_min), np.log10(cfg.O_max), 200)
    K = K_from_O(O, cfg.psi, cfg.c, cfg.alpha)
    F = Fstar(K, cfg.Fmax, cfg.K50)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    axes[0].plot(K, F, color="black", lw=2)
    axes[0].set(xscale="log", xlabel="competition ratio $K$",
                ylabel="evolved filter $F^*(K)$",
                title="Maternal filter rises with competition")
    axes[1].plot(O, F, color="#1F77B4", lw=2)
    axes[1].set(xscale="log", xlabel="ovule number $O$", ylabel="evolved filter $F^*$",
                title=f"Few ovules → stricter pistil  ($\\alpha={cfg.alpha}<1$)")
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


# =============================================================================
# Main
# =============================================================================
def main() -> None:
    ap = argparse.ArgumentParser(description="Layer 2 maternal filter + asymmetry.")
    ap.add_argument("--outdir", default="layer2_outputs")
    args = ap.parse_args()
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)

    print("[validate] checking the asymmetry core (sign theorem, linearisation, hump) ...")
    validate_core()

    cfg = Layer2Config()
    df = simulate_pairs(cfg)
    reg = validate_interaction(df)

    df.to_csv(out / "layer2_pairs.csv", index=False)
    reg.to_csv(out / "layer2_interaction_regression.csv", index=False)
    fig_hump(out / "fig_layer2_hump.png")
    fig_interaction(df, out / "fig_layer2_interaction.png")
    fig_filter(cfg, out / "fig_layer2_filter.png")

    # compact readout of the directional claim
    few = df[df["dlogO"] > 0]
    print(f"\n[claim] of {len(few)} pairs where i is the lower-ovule (higher-K) mother, "
          f"{(few['A_sobelchen'] > 0).mean()*100:.1f}% have A_ij > 0 (RI stronger into i).")
    print(f"[done] wrote pair table, regression, and 3 figures to: {out.resolve()}")


if __name__ == "__main__":
    main()
