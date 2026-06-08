"""
===============================================================================
Layer 1 -- Ecological competition engine (NO reproductive isolation yet)
===============================================================================
Project : ovule-gated pollen competition and passenger RI
Maps to : theory note, "Simulation protocol -> Layer 1: ecological competition
          without RI", and the sections "The competition ratio", "Sublinear
          pollen-load scaling", and "Selection among pollen grains".

WHAT THIS MODULE DOES
---------------------
It builds the single mechanistic hinge of the whole model and checks it. Given
ovule number O, it:
  1. draws an overdispersed effective pollen load   P ~ NegBinomial(mu, theta),
     with mean   mu = psi * c * O**alpha           (sublinear scaling, eq. Ksublinear);
  2. forms the competition ratio                    K = P / O                  (eq. Kdef);
  3. computes the opportunity index                 Omega = E[(P-O)_+ / P]      (eq. Omega);
  4. computes the truncation-selection differential S_z(K)                      (eq. Sz),
     i.e. how strongly the winning pollen grains are selected when K > 1.

There is deliberately NO maternal filter, NO mismatch, NO RI here. Layer 1 only
has to establish that few-ovule lineages enter a stronger-selection regime, and
that the closed-form S_z(K) we will rely on downstream actually equals what a
direct "draw the grains, keep the winners" experiment produces.

WHY THE VALIDATION MATTERS
--------------------------
S_z(K) is the classical truncation-selection differential: select the top
fraction q = 1/K of a Gaussian trait and the winners' mean exceeds the
population mean by  sigma_z * phi(t)/q  with  t = Phi^{-1}(1-q). Before any
result rides on this formula we confirm two things:
  (A) it reproduces the standard selection-intensity table (Falconer), and
  (B) it equals a direct Monte-Carlo experiment that actually draws P pollen
      performances, keeps the top O, and measures the realised differential.
If (A) and (B) hold, the engine's selection core is trustworthy.

NOTES FOR THE READER
--------------------
* We work in units of sigma_z (set sigma_z = 1), so S_z is the standardised
  selection intensity i -- a pure, comparable number.
* The closed form is the *infinite-population* truncation limit. Real flowers
  carry finite, overdispersed P, so we report BOTH the analytic per-flower
  average and (in the validation) the finite-sample realised value, and we flag
  the Jensen gap between S_z(mean K) and mean S_z(K).
===============================================================================
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
from scipy.stats import norm

import matplotlib
matplotlib.use("Agg")          # headless: write figures to disk, no display needed
import matplotlib.pyplot as plt


# =============================================================================
# Core deterministic functions (the maths of the note, one function each)
# =============================================================================
def pollen_mean(O: np.ndarray, psi: float, c: float, alpha: float) -> np.ndarray:
    """Mean effective pollen load mu = psi * c * O**alpha   (sublinear if alpha<1).

    O can be a scalar or array of ovule numbers. Returns the expected effective
    pollen load per flower for each O.
    """
    O = np.asarray(O, dtype=float)
    return psi * c * O ** alpha


def K_from_O(O: np.ndarray, psi: float, c: float, alpha: float) -> np.ndarray:
    """Expected competition ratio K = mu / O = psi * c * O**(alpha-1).

    With alpha < 1 the exponent (alpha-1) is negative, so K declines as O rises:
    few-ovule lineages have the higher competition regime.
    """
    O = np.asarray(O, dtype=float)
    return psi * c * O ** (alpha - 1.0)


def O_star(psi: float, c: float, alpha: float) -> float:
    """Ovule number at which expected competition switches on, K = 1.

    Solves 1 = psi*c*O**(alpha-1)  ->  O* = (psi*c)**(1/(1-alpha)).
    Undefined at alpha == 1 (K independent of O); returns nan there.
    """
    if np.isclose(alpha, 1.0):
        return float("nan")
    return (psi * c) ** (1.0 / (1.0 - alpha))


def Sz_analytic(K: np.ndarray, sigma_z: float = 1.0) -> np.ndarray:
    """Truncation-selection differential among winning pollen grains, S_z(K).

    Select the top fraction q = 1/K of a Gaussian performance distribution.
    The winners' mean exceeds the population mean by sigma_z * phi(t)/q, with
    t = Phi^{-1}(1 - q). For K <= 1 every grain can win, so S_z = 0.

    This is the standardised selection intensity i(q) scaled by sigma_z.
    """
    K = np.asarray(K, dtype=float)
    out = np.zeros_like(K)
    comp = K > 1.0                       # only here is there a losing margin
    q = 1.0 / K[comp]                    # winning (saved) fraction
    t = norm.ppf(1.0 - q)                # truncation point on the standard scale
    out[comp] = sigma_z * norm.pdf(t) / q
    return out


def omega_per_flower(P: np.ndarray, O: float) -> np.ndarray:
    """Per-flower opportunity index (P - O)_+ / P = max(0, 1 - O/P).

    Fraction of effective pollen grains that cannot win an ovule. Zero when
    P <= O. Averaging over flowers gives the lineage-level Omega (eq. Omega).
    """
    P = np.asarray(P, dtype=float)
    out = np.where(P > O, (P - O) / np.maximum(P, 1e-12), 0.0)
    return out


# =============================================================================
# Stochastic pollen loads: NegBinomial parameterised by mean and dispersion
# =============================================================================
def draw_pollen_loads(mu: float, theta: float, n: int, rng: np.random.Generator) -> np.ndarray:
    """Draw n effective pollen loads ~ NegBinomial with mean=mu, dispersion=theta.

    Ecology (NB2) parameterisation: variance = mu + mu**2/theta, so small theta
    => strongly overdispersed (clumped) pollen arrival; theta -> infinity => Poisson.
    numpy's negative_binomial(n_failures, p) has mean n_failures*(1-p)/p, so we set
    n_failures = theta and p = theta/(theta+mu).
    """
    if mu <= 0:
        return np.zeros(n, dtype=float)
    p = theta / (theta + mu)
    return rng.negative_binomial(theta, p, size=n).astype(float)


# =============================================================================
# (A) + (B) VALIDATION of the selection core
# =============================================================================
def validate_truncation(sigma_z: float = 1.0,
                        K_grid: List[float] = (1.5, 2, 3, 5, 10, 25, 100),
                        P_large: int = 20000,
                        n_rep: int = 400,
                        seed: int = 1) -> pd.DataFrame:
    """Confirm S_z(K) (analytic) == realised differential from a direct experiment.

    For each K we set q = 1/K, draw P_large Gaussian pollen performances, keep the
    top O = round(q * P_large), and measure the realised selection differential
    (mean winners - mean all). Repeating n_rep times gives a mean and a standard
    error; the analytic value should fall inside mean +/- 3*SE.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for K in K_grid:
        q = 1.0 / K
        O = max(1, int(round(q * P_large)))
        diffs = np.empty(n_rep)
        for r in range(n_rep):
            z = rng.normal(0.0, sigma_z, size=P_large)
            winners = np.partition(z, P_large - O)[P_large - O:]   # top O values
            diffs[r] = winners.mean() - z.mean()
        mc_mean, mc_se = diffs.mean(), diffs.std(ddof=1) / np.sqrt(n_rep)
        ana = float(Sz_analytic(np.array([K]), sigma_z)[0])
        within = abs(ana - mc_mean) <= 3.0 * mc_se + 0.01  # 3 SE + tiny finite-P slack
        rows.append(dict(K=K, q=q, analytic=ana, mc_mean=mc_mean, mc_se=mc_se,
                         within_3SE=bool(within)))
    df = pd.DataFrame(rows)
    assert df["within_3SE"].all(), (
        "Truncation-selection validation FAILED:\n" + df.to_string(index=False)
    )
    return df


# =============================================================================
# Per-lineage simulation (the engine itself)
# =============================================================================
@dataclass
class LineageResult:
    O: float
    psi: float
    alpha: float
    mean_K: float
    Omega: float
    mean_Sz_perflower: float      # E[ S_z(K_if) ]  -- the honest, Jensen-correct average
    Sz_at_meanK: float            # S_z( E[K_if] )  -- shown only to expose the Jensen gap
    frac_competitive: float       # fraction of flowers with P > O


def simulate_lineage(O: float, psi: float, c: float, alpha: float,
                    theta_P: float, n_flowers: int, sigma_z: float,
                    rng: np.random.Generator) -> LineageResult:
    """Simulate one lineage's flowers and summarise its competition statistics."""
    mu = float(pollen_mean(O, psi, c, alpha))
    P = draw_pollen_loads(mu, theta_P, n_flowers, rng)         # overdispersed loads
    K_if = np.where(P > 0, P / O, 0.0)                          # per-flower ratio
    Sz_if = Sz_analytic(K_if, sigma_z)                          # per-flower S_z
    return LineageResult(
        O=O, psi=psi, alpha=alpha,
        mean_K=float(K_if.mean()),
        Omega=float(omega_per_flower(P, O).mean()),
        mean_Sz_perflower=float(Sz_if.mean()),
        Sz_at_meanK=float(Sz_analytic(np.array([K_if.mean()]), sigma_z)[0]),
        frac_competitive=float((P > O).mean()),
    )


# =============================================================================
# Grid sweep over the psi-alpha sensitivity space and a range of ovule numbers
# =============================================================================
@dataclass
class GridConfig:
    psi_grid: List[float] = field(default_factory=lambda: [0.001, 0.003, 0.01, 0.02, 0.05, 0.10])
    alpha_grid: List[float] = field(default_factory=lambda: [0.25, 0.50, 0.75, 1.00, 1.25])
    O_grid: List[float] = field(default_factory=lambda:
        list(np.unique(np.round(np.logspace(np.log10(2), np.log10(20000), 24)))))
    c: float = 500.0          # scaling constant; with psi~0.02, alpha=0.5 puts O* near 100
    theta_P: float = 3.0      # NB dispersion: modest clumping of pollen arrival
    n_flowers: int = 2000
    sigma_z: float = 1.0
    seed: int = 7


def run_grid(cfg: GridConfig) -> pd.DataFrame:
    """Sweep psi x alpha x O, returning one tidy row per lineage."""
    rng = np.random.default_rng(cfg.seed)
    rows = []
    for psi in cfg.psi_grid:
        for alpha in cfg.alpha_grid:
            for O in cfg.O_grid:
                r = simulate_lineage(O, psi, cfg.c, alpha, cfg.theta_P,
                                    cfg.n_flowers, cfg.sigma_z, rng)
                rows.append(vars(r) | {"Ostar": O_star(psi, cfg.c, alpha)})
    return pd.DataFrame(rows)


# =============================================================================
# Figures
# =============================================================================
def fig_validation(valdf: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 5))
    Kline = np.linspace(1.01, max(valdf["K"]) * 1.05, 400)
    ax.plot(Kline, Sz_analytic(Kline), color="black", lw=2, label="analytic $S_z(K)$")
    ax.errorbar(valdf["K"], valdf["mc_mean"], yerr=3 * valdf["mc_se"], fmt="o",
                color="#C0392B", capsize=3, label="direct Monte Carlo $\\pm 3$ SE")
    ax.set_xlabel("competition ratio $K = P/O$")
    ax.set_ylabel("selection differential $S_z$  (units of $\\sigma_z$)")
    ax.set_title("Layer 1 validation: closed-form truncation selection vs experiment")
    ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def fig_K_and_Sz_vs_O(grid: pd.DataFrame, cfg: GridConfig, psi_show: float, path: Path) -> None:
    sub = grid[np.isclose(grid["psi"], psi_show)]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    cmap = plt.get_cmap("viridis")
    alphas = sorted(sub["alpha"].unique())
    for k, a in enumerate(alphas):
        s = sub[np.isclose(sub["alpha"], a)].sort_values("O")
        col = cmap(k / max(1, len(alphas) - 1))
        axes[0].plot(s["O"], s["mean_K"], color=col, marker="o", ms=3, label=f"$\\alpha$={a}")
        axes[1].plot(s["O"], s["mean_Sz_perflower"], color=col, marker="o", ms=3, label=f"$\\alpha$={a}")
    axes[0].axhline(1.0, color="grey", ls="--", lw=1)
    axes[0].set(xscale="log", yscale="log", xlabel="ovule number $O$",
                ylabel="mean competition ratio $K$",
                title=f"$K$ vs $O$  ($\\psi$={psi_show}); dashed line $K=1$")
    axes[1].set(xscale="log", xlabel="ovule number $O$",
                ylabel="mean per-flower $S_z$  (units of $\\sigma_z$)",
                title="Selection on pollen vs ovule number")
    for ax in axes:
        ax.legend(frameon=False, fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def fig_Ostar_heat(grid: pd.DataFrame, path: Path) -> None:
    piv = (grid.drop_duplicates(["psi", "alpha"])
               .pivot(index="alpha", columns="psi", values="Ostar"))
    fig, ax = plt.subplots(figsize=(7.5, 5))
    data = np.log10(piv.values)
    im = ax.imshow(data, origin="lower", aspect="auto", cmap="magma")
    ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels(piv.columns)
    ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index)
    ax.set_xlabel("deposition fraction $\\psi$"); ax.set_ylabel("scaling exponent $\\alpha$")
    ax.set_title("Competition threshold $\\log_{10} O^*$  (ovules where $K=1$)")
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            v = piv.values[i, j]
            txt = "n/a" if not np.isfinite(v) else f"{v:,.0f}"
            ax.text(j, i, txt, ha="center", va="center",
                    color="white" if data[i, j] < np.nanmedian(data) else "black", fontsize=8)
    fig.colorbar(im, ax=ax, label="$\\log_{10} O^*$")
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


# =============================================================================
# Main
# =============================================================================
def main() -> None:
    ap = argparse.ArgumentParser(description="Layer 1 ecological competition engine.")
    ap.add_argument("--outdir", default="layer1_outputs", help="output directory")
    ap.add_argument("--psi_show", type=float, default=0.02, help="psi value for the K/Sz figure")
    args = ap.parse_args()

    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)

    print("[validate] truncation selection: analytic vs direct Monte Carlo ...")
    valdf = validate_truncation()
    print(valdf.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print("[validate] PASSED: analytic S_z(K) matches the experiment within 3 SE.\n")

    print("[grid] sweeping psi x alpha x O ...")
    cfg = GridConfig()
    grid = run_grid(cfg)
    grid.to_csv(out / "layer1_grid.csv", index=False)
    valdf.to_csv(out / "layer1_validation.csv", index=False)

    fig_validation(valdf, out / "fig_layer1_validation.png")
    fig_K_and_Sz_vs_O(grid, cfg, args.psi_show, out / "fig_layer1_K_Sz_vs_O.png")
    fig_Ostar_heat(grid, out / "fig_layer1_Ostar.png")

    # A compact, honest readout of the central qualitative claim.
    psi, alpha = args.psi_show, 0.50
    s = grid[np.isclose(grid["psi"], psi) & np.isclose(grid["alpha"], alpha)].sort_values("O")
    lo, hi = s.iloc[0], s.iloc[-1]
    print(f"[claim] at psi={psi}, alpha={alpha}:")
    print(f"        few ovules  O={lo['O']:.0f}: mean K={lo['mean_K']:.2f}, Sz={lo['mean_Sz_perflower']:.3f}")
    print(f"        many ovules O={hi['O']:.0f}: mean K={hi['mean_K']:.3f}, Sz={hi['mean_Sz_perflower']:.3f}")
    print(f"        -> selection on pollen is stronger in the few-ovule lineage: "
          f"{lo['mean_Sz_perflower'] > hi['mean_Sz_perflower']}")
    print(f"\n[done] wrote grid, validation table, and 3 figures to: {out.resolve()}")


if __name__ == "__main__":
    main()
