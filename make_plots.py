from __future__ import annotations

import os
from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd


def _maybe_read(path: str) -> Optional[pd.DataFrame]:
    return pd.read_csv(path) if os.path.exists(path) else None


def _savefig(path: str) -> None:
    plt.tight_layout()
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()


def generate_all_plots(output_dir: str = "outputs") -> None:
    raw_dir = os.path.join(output_dir, "raw")
    fig_dir = os.path.join(output_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    asian_raw = _maybe_read(os.path.join(raw_dir, "curvature_raw.csv"))
    asian_exp = _maybe_read(os.path.join(raw_dir, "scaling_exponents.csv"))
    euro_raw = _maybe_read(os.path.join(raw_dir, "european_curvature_raw.csv"))
    euro_exp = _maybe_read(os.path.join(raw_dir, "european_scaling_exponents.csv"))

    if asian_raw is not None and not asian_raw.empty:
        fig = plt.figure(figsize=(9, 5.5))
        ax = fig.add_subplot(111)
        for H, grp in asian_raw.groupby("H"):
            grp = grp.sort_values("T", ascending=False)
            ax.plot(grp["T"], grp["abs_curv_mall"], marker="o", label=f"Mall H={H:.2f}")
            ax.plot(grp["T"], grp["abs_curv_fd"], marker="x", linestyle="--", label=f"FD H={H:.2f}")
            ax.plot(grp["T"], grp["abs_curv_theory"], linestyle=":")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Maturity T")
        ax.set_ylabel("Absolute curvature")
        ax.set_title("Asian curvature scaling")
        ax.legend(ncol=2, fontsize=8)
        _savefig(os.path.join(fig_dir, "scaling_loglog.png"))

        fig = plt.figure(figsize=(9, 5.5))
        ax = fig.add_subplot(111)
        for H, grp in asian_raw.groupby("H"):
            grp = grp.sort_values("T", ascending=False)
            denom = grp["T"] ** (2.0 * H + 1.0)
            ax.plot(grp["T"], grp["curv_mall"] / denom, marker="o", label=f"Mall H={H:.2f}")
            ax.plot(grp["T"], grp["curv_fd"] / denom, marker="x", linestyle="--", label=f"FD H={H:.2f}")
        ax.set_xscale("log")
        ax.set_xlabel("Maturity T")
        ax.set_ylabel(r"$\kappa_A / T^{2H+1}$")
        ax.set_title("Normalized Asian curvature")
        ax.legend(ncol=2, fontsize=8)
        _savefig(os.path.join(fig_dir, "normalized_curvature.png"))

        fig = plt.figure(figsize=(8, 5))
        ax = fig.add_subplot(111)
        if asian_exp is not None and not asian_exp.empty:
            ax.plot(asian_exp["H"], asian_exp["beta_error_mall"], marker="o", label="Malliavin")
            ax.plot(asian_exp["H"], asian_exp["beta_error_fd"], marker="x", label="FD")
        ax.set_xlabel("H")
        ax.set_ylabel("Absolute slope error")
        ax.set_title("Asian exponent errors against theory")
        ax.legend()
        _savefig(os.path.join(fig_dir, "exponent_errors.png"))

        fig = plt.figure(figsize=(8, 5))
        ax = fig.add_subplot(111)
        tmp = asian_raw.copy()
        tmp["rel_diff"] = (tmp["curv_fd"] - tmp["curv_mall"]).abs() / tmp[["abs_curv_mall", "abs_curv_fd"]].max(axis=1).clip(lower=1e-15)
        for H, grp in tmp.groupby("H"):
            grp = grp.sort_values("T", ascending=False)
            ax.plot(grp["T"], grp["rel_diff"], marker="o", label=f"H={H:.2f}")
        ax.set_xscale("log")
        ax.set_xlabel("Maturity T")
        ax.set_ylabel("Relative discrepancy")
        ax.set_title("Asian Malliavin and FD relative discrepancies")
        ax.legend(fontsize=8)
        _savefig(os.path.join(fig_dir, "cross_method_discrepancy.png"))

    if asian_exp is not None and euro_exp is not None and not asian_exp.empty and not euro_exp.empty:
        merged = asian_exp.merge(euro_exp, on="H", how="inner")
        fig = plt.figure(figsize=(8, 5))
        ax = fig.add_subplot(111)
        ax.plot(merged["H"], merged["beta_mall"], marker="o", label="Asian Malliavin")
        ax.plot(merged["H"], merged["beta_fd_euro"], marker="s", label="European FD")
        ax.plot(merged["H"], 2.0 * merged["H"] + 1.0, linestyle="--", label="Asian theory 2H+1")
        ax.plot(merged["H"], 2.0 * merged["H"], linestyle=":", label="European theory 2H")
        ax.set_xlabel("H")
        ax.set_ylabel("Fitted slope")
        ax.set_title("European and Asian scaling exponents")
        ax.legend(fontsize=8)
        _savefig(os.path.join(fig_dir, "european_and_asian.png"))


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description="Generate repository figures from CSV outputs.")
    parser.add_argument("--output-dir", default="outputs")
    args = parser.parse_args()
    generate_all_plots(args.output_dir)


if __name__ == "__main__":
    main()
