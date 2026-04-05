from __future__ import annotations

import argparse
import copy
import os

from asian_validation_common import ValidationConfig, run_european_control, run_full_validation, write_csv, write_json


def _parse_int_list(text: str):
    return [int(x.strip()) for x in text.split(',') if x.strip()]


def _parse_float_list(text: str):
    return [float(x.strip()) for x in text.split(',') if x.strip()]


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run robustness sweeps for Asian and European controls.")
    p.add_argument("--mode", choices=["smoke", "full"], default="smoke")
    p.add_argument("--output-dir", default="outputs")
    p.add_argument("--n-steps-grid", type=_parse_int_list, default=[128, 256])
    p.add_argument("--fd-bump-grid", type=_parse_float_list, default=[0.10, 0.20, 0.40])
    p.add_argument("--path-multipliers", type=_parse_float_list, default=[0.5, 1.0])
    return p


def main() -> None:
    args = build_parser().parse_args()
    base = ValidationConfig.smoke(output_dir=args.output_dir) if args.mode == "smoke" else ValidationConfig(output_dir=args.output_dir, mode="full")
    rows = []
    robustness_dir = os.path.join(args.output_dir, "summaries")
    os.makedirs(robustness_dir, exist_ok=True)

    for n_steps in args.n_steps_grid:
        for mult in args.path_multipliers:
            cfg = copy.deepcopy(base)
            cfg.output_dir = args.output_dir
            cfg.n_steps = int(n_steps)
            cfg.fd_bump_grid = list(args.fd_bump_grid)
            cfg.n_paths_mall = max(200, int(round(cfg.n_paths_mall * mult)))
            cfg.n_paths_fd = max(200, int(round(cfg.n_paths_fd * mult)))
            asian = run_full_validation(cfg)
            euro = run_european_control(cfg)
            rows.append({
                "n_steps": cfg.n_steps,
                "path_multiplier": mult,
                "n_paths_mall": cfg.n_paths_mall,
                "n_paths_fd": cfg.n_paths_fd,
                "asian_max_abs_beta_error_mall": asian["summary"]["max_abs_beta_error_mall"],
                "asian_max_abs_beta_error_fd": asian["summary"]["max_abs_beta_error_fd"],
                "asian_fd_blowup_fraction": asian["summary"]["fd_blowup_fraction"],
                "euro_max_abs_beta_error_fd": euro["summary"]["max_abs_beta_error_euro_fd"],
                "euro_fd_blowup_fraction": euro["summary"]["fd_blowup_fraction"],
            })

    write_csv(os.path.join(robustness_dir, "robustness_summary.csv"), rows)
    write_json(os.path.join(robustness_dir, "robustness_summary.json"), {"rows": rows})


if __name__ == "__main__":
    main()
