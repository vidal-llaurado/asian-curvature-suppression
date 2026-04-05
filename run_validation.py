from __future__ import annotations

import argparse
import os

from asian_validation_common import ValidationConfig, run_full_validation


def _parse_float_list(text: str):
    return [float(x.strip()) for x in text.split(',') if x.strip()]


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run Asian curvature Monte Carlo validation.")
    p.add_argument("--mode", choices=["smoke", "full"], default="full")
    p.add_argument("--output-dir", default="outputs")
    p.add_argument("--H-values", type=_parse_float_list)
    p.add_argument("--T-values", type=_parse_float_list)
    p.add_argument("--n-paths-mall", type=int)
    p.add_argument("--n-paths-fd", type=int)
    p.add_argument("--n-steps", type=int)
    p.add_argument("--n-seeds-fd", type=int)
    p.add_argument("--sigma0", type=float)
    p.add_argument("--nu", type=float)
    p.add_argument("--rho", type=float)
    p.add_argument("--fd-bump-frac", type=float)
    p.add_argument("--fd-bump-grid", type=_parse_float_list)
    p.add_argument("--skip-plots", action="store_true")
    return p


def build_config(args: argparse.Namespace) -> ValidationConfig:
    cfg = ValidationConfig.smoke(output_dir=args.output_dir) if args.mode == "smoke" else ValidationConfig(output_dir=args.output_dir, mode="full")
    if args.H_values is not None:
        cfg.H_values = args.H_values
    if args.T_values is not None:
        cfg.T_values = args.T_values
    if args.n_paths_mall is not None:
        cfg.n_paths_mall = args.n_paths_mall
    if args.n_paths_fd is not None:
        cfg.n_paths_fd = args.n_paths_fd
    if args.n_steps is not None:
        cfg.n_steps = args.n_steps
    if args.n_seeds_fd is not None:
        cfg.n_seeds_fd = args.n_seeds_fd
    if args.sigma0 is not None:
        cfg.sigma0 = args.sigma0
    if args.nu is not None:
        cfg.nu = args.nu
    if args.rho is not None:
        cfg.rho = args.rho
    if args.fd_bump_frac is not None:
        cfg.fd_bump_frac = args.fd_bump_frac
    if args.fd_bump_grid is not None:
        cfg.fd_bump_grid = args.fd_bump_grid
    return cfg


def main() -> None:
    args = build_parser().parse_args()
    cfg = build_config(args)
    run_full_validation(cfg)
    if not args.skip_plots:
        try:
            from make_plots import generate_all_plots
            generate_all_plots(cfg.output_dir)
        except Exception as exc:
            print(f"Plot generation failed: {exc}")


if __name__ == "__main__":
    main()
