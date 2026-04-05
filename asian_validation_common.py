from __future__ import annotations

import csv
import json
import math
import os
import time
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import brentq
from scipy.special import beta as beta_func
from scipy.stats import linregress, norm

try:
    from numba import njit
    NUMBA_AVAILABLE = True
except Exception:
    NUMBA_AVAILABLE = False

    def njit(*args, **kwargs):
        def wrapper(func):
            return func
        return wrapper


@dataclass
class ValidationConfig:
    H_values: List[float] = field(default_factory=lambda: [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40])
    T_values: List[float] = field(default_factory=lambda: [0.1, 0.05, 0.025, 0.0125, 0.00625])
    n_paths_mall: int = 5_000
    n_paths_fd: int = 2_500
    n_steps: int = 256
    n_seeds_fd: int = 5
    sigma0: float = 0.3
    nu: float = 0.5
    rho: float = -0.7
    fd_bump_frac: float = 0.10
    fd_bump_grid: List[float] = field(default_factory=lambda: [0.10, 0.20, 0.40])
    fd_smile_shifts: List[int] = field(default_factory=lambda: [-2, -1, 0, 1, 2])
    fd_local_poly_degree: int = 2
    fd_stability_rel_se_max: float = 0.75
    fd_stability_neighbor_logtol: float = 1.10
    antithetic: bool = True
    rng_seed_mall: int = 42
    fd_seed_base: int = 42
    use_numba: bool = True
    output_dir: str = "."
    mode: str = "full"

    @staticmethod
    def smoke(output_dir: str = ".") -> "ValidationConfig":
        return ValidationConfig(
            H_values=[0.10, 0.30],
            T_values=[0.1, 0.05, 0.025, 0.0125],
            n_paths_mall=1_000,
            n_paths_fd=800,
            n_steps=128,
            n_seeds_fd=3,
            output_dir=output_dir,
            mode="smoke",
        )


THEORETICAL_C_H = {
    0.05: None,
    0.10: 1.138350,
    0.15: None,
    0.20: 0.686152,
    0.25: None,
    0.30: 0.436706,
    0.35: None,
    0.40: 0.290291,
    0.50: 0.200000,
}


def compute_C_H_numerically(H: float, n_quad: int = 200) -> float:
    b1 = beta_func(H + 0.5, 3.0)
    b2 = beta_func(H + 0.5, H + 4.5)
    from numpy.polynomial.legendre import leggauss

    nodes, weights = leggauss(n_quad)
    xi_nodes = 0.5 * (nodes + 1.0)
    xi_weights = 0.5 * weights
    t_nodes = 0.5 * (nodes + 1.0)
    t_weights = 0.5 * weights

    K = 0.0
    for xi, wxi in zip(xi_nodes, xi_weights):
        one_minus_xi = 1.0 - xi
        outer = one_minus_xi ** (H + 3.5)
        for t, wt in zip(t_nodes, t_weights):
            integrand = outer * (t ** (H - 0.5)) * ((1.0 - t) ** 2) * ((xi + one_minus_xi * t) ** (H - 0.5))
            K += wxi * wt * integrand
    return float(b1 * b2 + 2.0 * K)


def get_C_H(H: float) -> float:
    if H in THEORETICAL_C_H and THEORETICAL_C_H[H] is not None:
        return float(THEORETICAL_C_H[H])
    return compute_C_H_numerically(H)


def theoretical_prefactor(H: float, sigma0: float, nu: float, rho: float) -> float:
    C = get_C_H(H)
    return float((3.0 * np.sqrt(6.0 * np.pi) * rho**2 * nu**2 * H * sigma0) / (2.0 * H + 2.5) * C)


def theoretical_curvature(H: float, T: float, sigma0: float, nu: float, rho: float) -> float:
    return float(theoretical_prefactor(H, sigma0, nu, rho) * T ** (2.0 * H + 1.0))


def build_cholesky_factor(H: float, n_steps: int, dt: float) -> np.ndarray:
    def fbm_cov(s: float, t: float, H_: float) -> float:
        return 0.5 * (abs(s) ** (2 * H_) + abs(t) ** (2 * H_) - abs(t - s) ** (2 * H_))

    incr_cov = np.zeros((n_steps, n_steps), dtype=np.float64)
    for i in range(n_steps):
        ti0, ti1 = i * dt, (i + 1) * dt
        for j in range(n_steps):
            tj0, tj1 = j * dt, (j + 1) * dt
            incr_cov[i, j] = (
                fbm_cov(ti1, tj1, H)
                - fbm_cov(ti1, tj0, H)
                - fbm_cov(ti0, tj1, H)
                + fbm_cov(ti0, tj0, H)
            )
    incr_cov += 1e-14 * np.eye(n_steps)
    return np.linalg.cholesky(incr_cov)


class DeterministicCache:
    def __init__(self) -> None:
        self._chol: Dict[Tuple[float, int, float], np.ndarray] = {}
        self._kernel_pack: Dict[Tuple[float, int, float], Dict[str, np.ndarray]] = {}

    def cholesky(self, H: float, n_steps: int, dt: float) -> np.ndarray:
        key = (float(H), int(n_steps), float(dt))
        if key not in self._chol:
            self._chol[key] = build_cholesky_factor(H, n_steps, dt)
        return self._chol[key]

    def kernel_pack(self, H: float, T: float, n_steps: int) -> Dict[str, np.ndarray]:
        dt = T / n_steps
        key = (float(H), int(n_steps), float(dt))
        if key not in self._kernel_pack:
            self._kernel_pack[key] = precompute_kernel_pack(H=H, T=T, n_steps=n_steps)
        return self._kernel_pack[key]


def _rough_kernel_from_diff(diff_idx: int, dt: float, H: float) -> float:
    if diff_idx == 0:
        return (dt ** (H - 0.5)) / (H + 0.5)
    lag = diff_idx * dt
    if lag < 0.5 * dt:
        return (dt ** (H - 0.5)) / (H + 0.5)
    return lag ** (H - 0.5)


def precompute_kernel_pack(H: float, T: float, n_steps: int) -> Dict[str, np.ndarray]:
    dt = T / n_steps
    idx = np.arange(n_steps, dtype=np.int64)
    t_mid = (idx + 0.5) * dt

    kdiff = np.empty(n_steps, dtype=np.float64)
    for d in range(n_steps):
        kdiff[d] = _rough_kernel_from_diff(d, dt, H)

    lower = np.zeros((n_steps, n_steps), dtype=np.float64)
    for u in range(n_steps):
        for r in range(u + 1):
            lower[u, r] = kdiff[u - r]

    phi_weight = (T - t_mid) / T
    variance_weight = ((T - t_mid) ** 2) / (T ** 2)
    active_mask = np.zeros(n_steps, dtype=bool)
    if n_steps >= 2:
        active_mask[: n_steps - 2] = True

    return {
        "dt": np.array(dt, dtype=np.float64),
        "t_mid": t_mid,
        "lower": lower,
        "phi_weight": phi_weight,
        "variance_weight": variance_weight,
        "active_mask": active_mask,
    }


def generate_fbm_and_z(L: np.ndarray, n_paths: int, rng: np.random.Generator) -> Tuple[np.ndarray, np.ndarray]:
    n_steps = L.shape[0]
    Z = rng.standard_normal((n_paths, n_steps))
    dBH = Z @ L.T
    return dBH, Z


def simulate_fbergomi(
    H: float,
    T: float,
    n_paths: int,
    n_steps: int,
    sigma0: float,
    nu: float,
    rho: float,
    *,
    L: np.ndarray,
    rng: np.random.Generator,
    antithetic: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    dt = T / n_steps
    t_grid = np.linspace(0.0, T, n_steps + 1)

    half = n_paths // 2 if antithetic else n_paths
    dBH, Z = generate_fbm_and_z(L, half, rng)
    if antithetic:
        dBH = np.vstack([dBH, -dBH])
        Z = np.vstack([Z, -Z])
        n_eff = 2 * half
    else:
        n_eff = half

    BH = np.zeros((n_eff, n_steps + 1), dtype=np.float64)
    BH[:, 1:] = np.cumsum(dBH, axis=1)
    V = sigma0**2 * np.exp(nu * np.sqrt(2.0 * H) * BH - 0.5 * nu**2 * t_grid ** (2.0 * H))
    sigma = np.sqrt(np.maximum(V, 1e-20))
    phi = sigma * ((T - t_grid) / T)
    return sigma, phi, t_grid, dBH, Z


def simulate_bachelier_paths_from_sigma(
    sigma: np.ndarray,
    Z: np.ndarray,
    dt: float,
    rho: float,
    rng: np.random.Generator,
) -> np.ndarray:
    dW_prime = Z * np.sqrt(dt)
    dW_perp = rng.standard_normal(Z.shape) * np.sqrt(dt)
    dW = rho * dW_prime + np.sqrt(max(1.0 - rho**2, 0.0)) * dW_perp
    S = np.zeros((sigma.shape[0], sigma.shape[1]), dtype=np.float64)
    for i in range(sigma.shape[1] - 1):
        S[:, i + 1] = S[:, i] + sigma[:, i] * dW[:, i]
    return S


@njit(cache=True)
def _compute_lambda2_numba(
    sigma: np.ndarray,
    phi: np.ndarray,
    lower: np.ndarray,
    variance_weight: np.ndarray,
    phi_weight: np.ndarray,
    dt: float,
    c1: float,
    rho: float,
) -> np.ndarray:
    n_paths, n_cols = sigma.shape
    n_steps = n_cols
    lambda2 = np.zeros((n_paths, n_steps), dtype=np.float64)
    inner_A = np.zeros((n_paths, n_steps), dtype=np.float64)

    for p in range(n_paths):
        for r in range(n_steps):
            acc = 0.0
            for u in range(r, n_steps):
                acc += variance_weight[u] * sigma[p, u] * sigma[p, u] * lower[u, r]
            inner_A[p, r] = 2.0 * rho * c1 * dt * acc

    for p in range(n_paths):
        for s in range(n_steps - 2):
            integral_r = 0.0
            for r in range(s, n_steps - 1):
                base_B = 0.0
                for u in range(r, n_steps):
                    base_B += variance_weight[u] * sigma[p, u] * sigma[p, u] * lower[u, r] * lower[u, s]
                inner_B = 4.0 * rho * rho * c1 * c1 * dt * base_B
                D_sW_phi_r = rho * phi_weight[r] * sigma[p, r] * c1 * lower[r, s]
                comp_A = D_sW_phi_r * inner_A[p, r]
                comp_B = phi[p, r] * inner_B
                integral_r += (comp_A + comp_B) * dt
            lambda2[p, s] = phi[p, s] * integral_r
    return lambda2


def compute_lambda2_pathwise(
    sigma: np.ndarray,
    phi: np.ndarray,
    *,
    H: float,
    nu: float,
    rho: float,
    kernel_pack: Dict[str, np.ndarray],
    use_numba: bool = True,
) -> np.ndarray:
    dt = float(kernel_pack["dt"])
    lower = kernel_pack["lower"]
    variance_weight = kernel_pack["variance_weight"]
    phi_weight = kernel_pack["phi_weight"]
    c1 = nu * np.sqrt(2.0 * H) / 2.0

    sigma_mid = np.ascontiguousarray(sigma[:, :-1])
    phi_mid = np.ascontiguousarray(phi[:, :-1])

    if use_numba and NUMBA_AVAILABLE:
        return _compute_lambda2_numba(sigma_mid, phi_mid, lower, variance_weight, phi_weight, dt, c1, rho)

    weighted_sigma2 = (sigma_mid ** 2) * variance_weight[None, :]
    inner_A = 2.0 * rho * c1 * dt * (weighted_sigma2 @ lower)

    n_paths, n_steps = sigma_mid.shape
    lambda2 = np.zeros((n_paths, n_steps), dtype=np.float64)
    for s in range(n_steps - 2):
        weighted_by_s = weighted_sigma2 * lower[:, s][None, :]
        base_B = weighted_by_s @ lower
        inner_B = 4.0 * rho**2 * c1**2 * dt * base_B
        D_sW_phi_r = rho * phi_weight * sigma_mid * c1 * lower[:, s][None, :]
        comp = D_sW_phi_r * inner_A + phi_mid * inner_B
        lambda2[:, s] = phi_mid[:, s] * np.sum(comp[:, s:n_steps - 1], axis=1) * dt
    return lambda2


def compute_curvature_correction(lambda2: np.ndarray, T: float, sigma0: float, kernel_pack: Dict[str, np.ndarray]) -> Tuple[float, float, np.ndarray]:
    dt = float(kernel_pack["dt"])
    t_mid = kernel_pack["t_mid"]
    active_mask = kernel_pack["active_mask"]
    tau = T - t_mid
    v_s = sigma0 * tau / (T * np.sqrt(3.0))
    G_s = 1.0 / (v_s**3 * np.sqrt(tau))
    G_s[~active_mask] = 0.0

    integral = np.sum(lambda2 * G_s[None, :], axis=1) * dt
    vega_A = np.sqrt(T / (2.0 * np.pi))
    samples = integral / vega_A
    mean_corr = float(np.mean(samples))
    se_corr = float(np.std(samples, ddof=1) / np.sqrt(samples.shape[0])) if samples.shape[0] > 1 else 0.0
    return mean_corr, se_corr, samples


def bachelier_price(x: float, k: float, sigma: float, T: float) -> float:
    if sigma * np.sqrt(T) < 1e-15:
        return float(max(x - k, 0.0))
    d = (x - k) / (sigma * np.sqrt(T))
    return float(sigma * np.sqrt(T) * (d * norm.cdf(d) + norm.pdf(d)))


def bachelier_iv_inversion(price: float, x: float, k: float, T: float, iv_guess: float) -> float:
    intrinsic = max(x - k, 0.0)
    if price <= intrinsic + 1e-14:
        return 1e-8
    target = max(price, intrinsic + 1e-14)

    def obj(sig: float) -> float:
        return bachelier_price(x, k, sig, T) - target

    moneyness_term = abs(x - k) / max(np.sqrt(T), 1e-12)
    extrinsic_term = max(target - intrinsic, 1e-14) * np.sqrt(2.0 * np.pi / max(T, 1e-16))
    hi = max(10.0 * max(iv_guess, 1e-4), moneyness_term + extrinsic_term + 1.0)
    lo = 1e-8

    f_lo = obj(lo)
    f_hi = obj(hi)
    if f_lo > 0.0:
        return lo
    if f_hi < 0.0:
        for _ in range(8):
            hi *= 2.0
            f_hi = obj(hi)
            if f_hi >= 0.0:
                break
        else:
            return hi
    try:
        return float(brentq(obj, lo, hi, xtol=1e-12, maxiter=200))
    except Exception:
        return float(max(iv_guess, lo))


def _parabolic_fit(shifts: Sequence[int], values_by_shift: Dict[int, float], dk: float) -> Tuple[np.ndarray, Dict[str, float]]:
    shifts_arr = np.asarray(list(shifts), dtype=np.float64)
    values = np.asarray([values_by_shift[int(s)] for s in shifts_arr], dtype=np.float64)
    X = np.column_stack([np.ones_like(shifts_arr), shifts_arr, shifts_arr**2])
    bandwidth = max(float(np.max(np.abs(shifts_arr))), 1.0) + 1e-12
    u = np.abs(shifts_arr) / bandwidth
    local_w = (1.0 - u**3) ** 3
    local_w = np.where(np.isfinite(local_w) & (local_w > 1e-10), local_w, 1e-10)
    XtWX = X.T @ (local_w[:, None] * X)
    XtWy = X.T @ (local_w * values)
    try:
        coef = np.linalg.solve(XtWX, XtWy)
    except np.linalg.LinAlgError:
        coef = np.linalg.lstsq(XtWX + 1e-12 * np.eye(3), XtWy, rcond=None)[0]
    fitted = X @ coef
    resid = values - fitted
    return coef, {
        "fit_rmse": float(np.sqrt(np.mean(resid**2))),
        "fit_cond": float(np.linalg.cond(XtWX)),
        "atm_fit": float(coef[0]),
        "linear_coef_scaled": float(coef[1] / dk),
        "quadratic_coef_scaled": float(coef[2] / (dk**2)),
    }


def _aggregate_seed_curvatures(seed_records: List[Dict[str, Any]]) -> Dict[str, float]:
    if not seed_records:
        return {"curvature": np.nan, "se": np.nan, "rmse_mean": np.nan, "cond_median": np.nan}
    curvs = np.array([row["curvature"] for row in seed_records], dtype=np.float64)
    rmses = np.array([row.get("iv_fit_rmse", row.get("fit_rmse", np.nan)) for row in seed_records], dtype=np.float64)
    conds = np.array([row.get("iv_fit_cond", row.get("fit_cond", np.nan)) for row in seed_records], dtype=np.float64)
    weights = 1.0 / np.maximum(rmses, 1e-8) ** 2
    weights = np.where(np.isfinite(weights), weights, 0.0)
    if np.sum(weights) <= 0:
        weights = np.ones_like(curvs)
    weights = weights / np.sum(weights)
    mean_curv = float(np.sum(weights * curvs))
    centered = curvs - mean_curv
    eff_n = float(1.0 / np.sum(weights**2))
    var = float(np.sum(weights * centered**2))
    se = float(np.sqrt(max(var, 0.0) / max(eff_n, 1.0)))
    return {
        "curvature": mean_curv,
        "se": se,
        "rmse_mean": float(np.nanmean(rmses)),
        "cond_median": float(np.nanmedian(conds)),
        "effective_n": eff_n,
    }


def _curvature_from_price_smile(
    prices_by_shift: Dict[int, float],
    strike_shifts: Sequence[int],
    dk: float,
    T: float,
    iv_guess: float,
) -> Tuple[float, Dict[str, float], Dict[int, float]]:
    price_coef, price_diag = _parabolic_fit(strike_shifts, prices_by_shift, dk)
    smoothed_prices: Dict[int, float] = {}
    ivs_by_shift: Dict[int, float] = {}
    for shift in strike_shifts:
        z = float(shift)
        smoothed_price = float(price_coef[0] + price_coef[1] * z + price_coef[2] * z * z)
        k = shift * dk
        intrinsic = max(-k, 0.0)
        smoothed_price = max(smoothed_price, intrinsic + 1e-12)
        smoothed_prices[int(shift)] = smoothed_price
        ivs_by_shift[int(shift)] = bachelier_iv_inversion(smoothed_price, 0.0, k, T, iv_guess)

    iv_coef, iv_diag = _parabolic_fit(strike_shifts, ivs_by_shift, dk)
    curvature = float(2.0 * iv_coef[2] / (dk**2))
    diag = {
        "price_fit_rmse": price_diag["fit_rmse"],
        "price_fit_cond": price_diag["fit_cond"],
        "price_atm_fit": price_diag["atm_fit"],
        "price_quadratic_coef_scaled": price_diag["quadratic_coef_scaled"],
        "iv_fit_rmse": iv_diag["fit_rmse"],
        "iv_fit_cond": iv_diag["fit_cond"],
        "iv_atm_fit": iv_diag["atm_fit"],
        "iv_quadratic_coef_scaled": iv_diag["quadratic_coef_scaled"],
    }
    return curvature, diag, ivs_by_shift


def _asian_average(S: np.ndarray, dt: float, T: float) -> np.ndarray:
    n_steps = S.shape[1] - 1
    weights_trap = np.ones(n_steps + 1, dtype=np.float64)
    weights_trap[0] = 0.5
    weights_trap[-1] = 0.5
    return np.sum(S * weights_trap[None, :], axis=1) * dt / T


def finite_difference_curvature_asian(
    H: float,
    T: float,
    n_paths: int,
    n_steps: int,
    sigma0: float,
    nu: float,
    rho: float,
    *,
    n_seeds: int,
    fd_bump_frac: float,
    L: np.ndarray,
    antithetic: bool,
    fd_seed_base: int,
    strike_shifts: Sequence[int] = (-2, -1, 0, 1, 2),
) -> Tuple[float, float, List[Dict[str, float]]]:
    dt = T / n_steps
    v_asian_0 = sigma0 / np.sqrt(3.0)
    dk = v_asian_0 * np.sqrt(T) * fd_bump_frac
    strike_shifts = tuple(int(x) for x in strike_shifts)
    per_seed: List[Dict[str, float]] = []

    for seed in range(n_seeds):
        rng = np.random.default_rng(1000 * seed + fd_seed_base)
        sigma, _, _, _, Z = simulate_fbergomi(H, T, n_paths, n_steps, sigma0, nu, rho, L=L, rng=rng, antithetic=antithetic)
        S = simulate_bachelier_paths_from_sigma(sigma, Z, dt, rho, rng)
        A_T = _asian_average(S, dt, T)

        prices_by_shift: Dict[int, float] = {}
        record: Dict[str, float] = {"seed": float(seed), "dk": float(dk)}
        for shift in strike_shifts:
            k = shift * dk
            payoff = np.maximum(A_T - k, 0.0)
            price = float(np.mean(payoff))
            price_se = float(np.std(payoff, ddof=1) / np.sqrt(len(payoff))) if len(payoff) > 1 else 0.0
            prices_by_shift[int(shift)] = price
            key = f"m{abs(shift)}" if shift < 0 else (f"p{shift}" if shift > 0 else "0")
            record[f"price_{key}"] = price
            record[f"price_se_{key}"] = price_se

        curvature, diag, ivs_by_shift = _curvature_from_price_smile(prices_by_shift, strike_shifts, dk, T, v_asian_0)
        record["curvature"] = curvature
        record.update(diag)
        for shift in strike_shifts:
            key = f"m{abs(shift)}" if shift < 0 else (f"p{shift}" if shift > 0 else "0")
            record[f"iv_{key}"] = ivs_by_shift[int(shift)]
        per_seed.append(record)

    agg = _aggregate_seed_curvatures(per_seed)
    return float(agg["curvature"]), float(agg["se"]), per_seed


def finite_difference_curvature_european(
    H: float,
    T: float,
    n_paths: int,
    n_steps: int,
    sigma0: float,
    nu: float,
    rho: float,
    *,
    n_seeds: int,
    fd_bump_frac: float,
    L: np.ndarray,
    antithetic: bool,
    fd_seed_base: int,
    strike_shifts: Sequence[int] = (-2, -1, 0, 1, 2),
) -> Tuple[float, float, List[Dict[str, float]]]:
    dt = T / n_steps
    iv_guess = sigma0
    dk = iv_guess * np.sqrt(T) * fd_bump_frac
    strike_shifts = tuple(int(x) for x in strike_shifts)
    per_seed: List[Dict[str, float]] = []

    for seed in range(n_seeds):
        rng = np.random.default_rng(1000 * seed + fd_seed_base)
        sigma, _, _, _, Z = simulate_fbergomi(H, T, n_paths, n_steps, sigma0, nu, rho, L=L, rng=rng, antithetic=antithetic)
        S = simulate_bachelier_paths_from_sigma(sigma, Z, dt, rho, rng)
        S_T = S[:, -1]

        prices_by_shift: Dict[int, float] = {}
        record: Dict[str, float] = {"seed": float(seed), "dk": float(dk)}
        for shift in strike_shifts:
            k = shift * dk
            payoff = np.maximum(S_T - k, 0.0)
            price = float(np.mean(payoff))
            price_se = float(np.std(payoff, ddof=1) / np.sqrt(len(payoff))) if len(payoff) > 1 else 0.0
            prices_by_shift[int(shift)] = price
            key = f"m{abs(shift)}" if shift < 0 else (f"p{shift}" if shift > 0 else "0")
            record[f"price_{key}"] = price
            record[f"price_se_{key}"] = price_se

        curvature, diag, ivs_by_shift = _curvature_from_price_smile(prices_by_shift, strike_shifts, dk, T, iv_guess)
        record["curvature"] = curvature
        record.update(diag)
        for shift in strike_shifts:
            key = f"m{abs(shift)}" if shift < 0 else (f"p{shift}" if shift > 0 else "0")
            record[f"iv_{key}"] = ivs_by_shift[int(shift)]
        per_seed.append(record)

    agg = _aggregate_seed_curvatures(per_seed)
    return float(agg["curvature"]), float(agg["se"]), per_seed


def finite_difference_curvature_sweep(
    kind: str,
    H: float,
    T: float,
    n_paths: int,
    n_steps: int,
    sigma0: float,
    nu: float,
    rho: float,
    *,
    n_seeds: int,
    fd_bump_grid: Sequence[float],
    fd_bump_frac_default: float,
    L: np.ndarray,
    antithetic: bool,
    fd_seed_base: int,
    rel_se_max: float,
    neighbor_logtol: float,
    strike_shifts: Sequence[int] = (-2, -1, 0, 1, 2),
) -> Dict[str, Any]:
    estimator = finite_difference_curvature_asian if kind == "asian" else finite_difference_curvature_european
    bump_grid = sorted(float(x) for x in fd_bump_grid)
    bump_results: List[Dict[str, Any]] = []

    for bump in bump_grid:
        mean_curv, se_curv, per_seed = estimator(
            H, T, n_paths, n_steps, sigma0, nu, rho,
            n_seeds=n_seeds, fd_bump_frac=bump, L=L, antithetic=antithetic,
            fd_seed_base=fd_seed_base, strike_shifts=strike_shifts
        )
        abs_mean = abs(mean_curv)
        rel_se = se_curv / max(abs_mean, 1e-30)
        bump_results.append({
            "fd_bump_frac": bump,
            "mean_curv": mean_curv,
            "se_curv": se_curv,
            "rel_se": rel_se,
            "fit_rmse_mean": float(np.nanmean([row.get("iv_fit_rmse", np.nan) for row in per_seed])),
            "fit_cond_median": float(np.nanmedian([row.get("iv_fit_cond", np.nan) for row in per_seed])),
            "per_seed": per_seed,
        })

    for i, res in enumerate(bump_results):
        neighbors = []
        if i > 0:
            neighbors.append(bump_results[i - 1])
        if i + 1 < len(bump_results):
            neighbors.append(bump_results[i + 1])
        if neighbors and abs(res["mean_curv"]) > 1e-30:
            log_self = math.log(abs(res["mean_curv"]))
            res["neighbor_log_gap"] = float(max(abs(log_self - math.log(max(abs(n["mean_curv"]), 1e-30))) for n in neighbors))
        else:
            res["neighbor_log_gap"] = np.nan
        res["stable"] = bool(res["rel_se"] <= rel_se_max and (not np.isfinite(res["neighbor_log_gap"]) or res["neighbor_log_gap"] <= neighbor_logtol))

    stable = [r for r in bump_results if r["stable"]]
    if stable:
        selected = min(stable, key=lambda r: (r["rel_se"], abs(r["fd_bump_frac"] - fd_bump_frac_default)))
        fd_status = "stable"
    else:
        selected = min(bump_results, key=lambda r: (r["rel_se"], abs(r["fd_bump_frac"] - fd_bump_frac_default)))
        fd_status = "unstable_fallback"

    smallest = bump_results[0]
    largest = bump_results[-1]
    blowup_ratio = float(abs(smallest["mean_curv"]) / max(abs(largest["mean_curv"]), 1e-30))
    blowup_flag = bool(smallest["rel_se"] > rel_se_max or blowup_ratio > 5.0 or fd_status != "stable")
    return {
        "mean_curv": selected["mean_curv"],
        "se_curv": selected["se_curv"],
        "selected_bump_frac": selected["fd_bump_frac"],
        "fd_status": fd_status,
        "fd_blowup_flag": blowup_flag,
        "fd_blowup_ratio_smallest_to_largest": blowup_ratio,
        "fd_rel_se_selected": selected["rel_se"],
        "fd_neighbor_log_gap_selected": selected["neighbor_log_gap"],
        "bump_results": bump_results,
        "per_seed_selected": selected["per_seed"],
    }


def weighted_loglog_fit(T_values: Sequence[float], y_abs_values: Sequence[float], y_se_values: Optional[Sequence[float]] = None) -> Dict[str, float]:
    T_arr = np.asarray(T_values, dtype=np.float64)
    y_arr = np.asarray(y_abs_values, dtype=np.float64)
    valid = np.isfinite(T_arr) & np.isfinite(y_arr) & (T_arr > 0.0) & (y_arr > 0.0)
    if valid.sum() < 3:
        return {"beta": np.nan, "beta_se": np.nan, "intercept": np.nan, "r2": np.nan, "n_points": int(valid.sum())}
    x = np.log(T_arr[valid])
    y = np.log(y_arr[valid])
    if y_se_values is None:
        fit = linregress(x, y)
        return {"beta": float(fit.slope), "beta_se": float(fit.stderr), "intercept": float(fit.intercept), "r2": float(fit.rvalue**2), "n_points": int(valid.sum())}

    se_arr = np.asarray(y_se_values, dtype=np.float64)[valid]
    rel = np.maximum(se_arr / np.maximum(y_arr[valid], 1e-30), 1e-8)
    w = 1.0 / (rel**2)
    X = np.column_stack([np.ones_like(x), x])
    W = np.diag(w)
    XtWX = X.T @ W @ X
    XtWy = X.T @ W @ y
    beta_hat = np.linalg.solve(XtWX, XtWy)
    resid = y - X @ beta_hat
    y_bar = np.average(y, weights=w)
    ss_res = float(np.sum(w * resid**2))
    ss_tot = float(np.sum(w * (y - y_bar) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    dof = max(len(x) - 2, 1)
    sigma2 = ss_res / dof
    cov = sigma2 * np.linalg.inv(XtWX)
    return {"beta": float(beta_hat[1]), "beta_se": float(np.sqrt(max(cov[1, 1], 0.0))), "intercept": float(beta_hat[0]), "r2": float(r2), "n_points": int(valid.sum())}


def summarize_validation(config: ValidationConfig, raw_rows: List[Dict[str, Any]], exponent_rows: List[Dict[str, Any]], prefactor_rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    beta_errors_m, beta_errors_f, rel_discrepancies, pref_ratio_m, pref_ratio_f = [], [], [], [], []
    fd_blowups, fd_selected_rel_se = [], []
    for row in exponent_rows:
        beta_theory = row["beta_theory_asian"]
        if np.isfinite(row["beta_mall"]):
            beta_errors_m.append(abs(row["beta_mall"] - beta_theory))
        if np.isfinite(row["beta_fd"]):
            beta_errors_f.append(abs(row["beta_fd"] - beta_theory))
    for row in raw_rows:
        denom = max(abs(row["curv_mall"]), abs(row["curv_fd"]), 1e-15)
        rel_discrepancies.append(abs(row["curv_mall"] - row["curv_fd"]) / denom)
        fd_blowups.append(bool(row.get("fd_blowup_flag", False)))
        if np.isfinite(row.get("fd_rel_se_selected", np.nan)):
            fd_selected_rel_se.append(float(row["fd_rel_se_selected"]))
    for row in prefactor_rows:
        if np.isfinite(row["pref_ratio_mall"]):
            pref_ratio_m.append(row["pref_ratio_mall"])
        if np.isfinite(row["pref_ratio_fd"]):
            pref_ratio_f.append(row["pref_ratio_fd"])

    summary = {
        "config": asdict(config),
        "max_abs_beta_error_mall": float(np.max(beta_errors_m)) if beta_errors_m else np.nan,
        "max_abs_beta_error_fd": float(np.max(beta_errors_f)) if beta_errors_f else np.nan,
        "mean_cross_method_rel_discrepancy": float(np.mean(rel_discrepancies)) if rel_discrepancies else np.nan,
        "max_cross_method_rel_discrepancy": float(np.max(rel_discrepancies)) if rel_discrepancies else np.nan,
        "median_pref_ratio_mall": float(np.median(pref_ratio_m)) if pref_ratio_m else np.nan,
        "median_pref_ratio_fd": float(np.median(pref_ratio_f)) if pref_ratio_f else np.nan,
        "fd_blowup_fraction": float(np.mean(fd_blowups)) if fd_blowups else np.nan,
        "fd_selected_rel_se_median": float(np.median(fd_selected_rel_se)) if fd_selected_rel_se else np.nan,
        "n_raw_rows": len(raw_rows),
        "n_exponent_rows": len(exponent_rows),
        "n_prefactor_rows": len(prefactor_rows),
    }
    summary["checks"] = {
        "mall_beta_error_le_0p20": bool(summary["max_abs_beta_error_mall"] <= 0.20) if np.isfinite(summary["max_abs_beta_error_mall"]) else False,
        "fd_beta_error_le_0p35": bool(summary["max_abs_beta_error_fd"] <= 0.35) if np.isfinite(summary["max_abs_beta_error_fd"]) else False,
        "mean_cross_method_rel_discrepancy_le_0p50": bool(summary["mean_cross_method_rel_discrepancy"] <= 0.50) if np.isfinite(summary["mean_cross_method_rel_discrepancy"]) else False,
        "fd_blowup_fraction_le_0p25": bool(summary["fd_blowup_fraction"] <= 0.25) if np.isfinite(summary["fd_blowup_fraction"]) else False,
    }
    summary["all_primary_checks_pass"] = bool(all(summary["checks"].values()))
    return summary


def summarize_european_control(config: ValidationConfig, raw_rows: List[Dict[str, Any]], exponent_rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    beta_errors = []
    blowups = []
    for row in exponent_rows:
        if np.isfinite(row["beta_fd_euro"]):
            beta_errors.append(abs(row["beta_fd_euro"] - row["beta_theory_euro"]))
    for row in raw_rows:
        blowups.append(bool(row.get("fd_blowup_flag", False)))
    summary = {
        "config": asdict(config),
        "max_abs_beta_error_euro_fd": float(np.max(beta_errors)) if beta_errors else np.nan,
        "median_abs_beta_error_euro_fd": float(np.median(beta_errors)) if beta_errors else np.nan,
        "fd_blowup_fraction": float(np.mean(blowups)) if blowups else np.nan,
        "n_raw_rows": len(raw_rows),
        "n_exponent_rows": len(exponent_rows),
    }
    summary["checks"] = {
        "euro_beta_error_le_0p35": bool(summary["max_abs_beta_error_euro_fd"] <= 0.35) if np.isfinite(summary["max_abs_beta_error_euro_fd"]) else False,
        "fd_blowup_fraction_le_0p25": bool(summary["fd_blowup_fraction"] <= 0.25) if np.isfinite(summary["fd_blowup_fraction"]) else False,
    }
    summary["all_primary_checks_pass"] = bool(all(summary["checks"].values()))
    return summary


def ensure_output_dirs(output_dir: str) -> None:
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, "raw"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "summaries"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "figures"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "diagnostics"), exist_ok=True)


def write_csv(path: str, rows: List[Dict[str, Any]], fieldnames: Optional[List[str]] = None) -> None:
    if fieldnames is None and rows:
        fieldnames = list(rows[0].keys())
    if fieldnames is None:
        return
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_json(path: str, payload: Dict[str, Any]) -> None:
    with open(path, "w") as f:
        json.dump(payload, f, indent=2, sort_keys=True)


def run_full_validation(config: ValidationConfig) -> Dict[str, Any]:
    ensure_output_dirs(config.output_dir)
    cache = DeterministicCache()
    raw_rows: List[Dict[str, Any]] = []
    exponent_rows: List[Dict[str, Any]] = []
    prefactor_rows: List[Dict[str, Any]] = []
    detail_rows: List[Dict[str, Any]] = []
    all_results: Dict[float, Dict[str, Any]] = {}

    print("=" * 100)
    print("Asian Curvature Suppression — Monte Carlo Validation")
    print("=" * 100)
    print(f"Mode            : {config.mode}")
    print(f"Malliavin paths : {config.n_paths_mall} (antithetic={config.antithetic})")
    print(f"FD paths        : {config.n_paths_fd} (antithetic={config.antithetic}) x {config.n_seeds_fd} seeds")
    print(f"Time steps      : {config.n_steps}")
    print(f"sigma0={config.sigma0}, nu={config.nu}, rho={config.rho}, fd_bump_frac={config.fd_bump_frac}")
    print(f"Numba enabled   : {config.use_numba and NUMBA_AVAILABLE}")
    print(f"Output dir      : {os.path.abspath(config.output_dir)}")
    print()

    for H in config.H_values:
        print(f"─── H = {H:.2f} ───")
        print(f"  {'T':>10} | {'Mall curv':>14} | {'SE_M':>10} | {'FD curv':>14} | {'SE_FD':>10} | {'Theory':>14} | {'t_M(s)':>8} | {'t_FD(s)':>8}")
        print("  " + "-" * 108)
        mall_vals, mall_ses, fd_vals, fd_ses = [], [], [], []

        for T in config.T_values:
            dt = T / config.n_steps
            L = cache.cholesky(H, config.n_steps, dt)
            kernel_pack = cache.kernel_pack(H, T, config.n_steps)

            t0 = time.perf_counter()
            rng_mall = np.random.default_rng(config.rng_seed_mall)
            sigma, phi, _, _, _ = simulate_fbergomi(H, T, config.n_paths_mall, config.n_steps, config.sigma0, config.nu, config.rho, L=L, rng=rng_mall, antithetic=config.antithetic)
            lambda2 = compute_lambda2_pathwise(sigma, phi, H=H, nu=config.nu, rho=config.rho, kernel_pack=kernel_pack, use_numba=config.use_numba)
            mall_mean, mall_se, _ = compute_curvature_correction(lambda2, T, config.sigma0, kernel_pack)
            time_mall = time.perf_counter() - t0

            t0 = time.perf_counter()
            fd_out = finite_difference_curvature_sweep(
                "asian", H, T, config.n_paths_fd, config.n_steps, config.sigma0, config.nu, config.rho,
                n_seeds=config.n_seeds_fd, fd_bump_grid=config.fd_bump_grid, fd_bump_frac_default=config.fd_bump_frac,
                L=L, antithetic=config.antithetic, fd_seed_base=config.fd_seed_base,
                rel_se_max=config.fd_stability_rel_se_max, neighbor_logtol=config.fd_stability_neighbor_logtol,
                strike_shifts=config.fd_smile_shifts,
            )
            fd_mean, fd_se, per_seed = fd_out["mean_curv"], fd_out["se_curv"], fd_out["per_seed_selected"]
            time_fd = time.perf_counter() - t0
            curv_th = theoretical_curvature(H, T, config.sigma0, config.nu, config.rho)

            abs_mall, abs_fd, abs_th = abs(mall_mean), abs(fd_mean), abs(curv_th)
            mall_vals.append(abs_mall)
            mall_ses.append(mall_se)
            fd_vals.append(abs_fd)
            fd_ses.append(fd_se)

            print(f"  {T:10.5f} | {mall_mean:14.6e} | {mall_se:10.2e} | {fd_mean:14.6e} | {fd_se:10.2e} | {curv_th:14.6e} | {time_mall:8.2f} | {time_fd:8.2f}")
            raw_rows.append({
                "H": H, "T": T, "logT": math.log(T),
                "curv_mall": mall_mean, "se_mall": mall_se, "abs_curv_mall": abs_mall,
                "log_abs_curv_mall": math.log(abs_mall) if abs_mall > 0 else np.nan,
                "curv_fd": fd_mean, "se_fd": fd_se, "abs_curv_fd": abs_fd,
                "log_abs_curv_fd": math.log(abs_fd) if abs_fd > 0 else np.nan,
                "curv_theory": curv_th, "abs_curv_theory": abs_th,
                "log_abs_curv_theory": math.log(abs_th) if abs_th > 0 else np.nan,
                "time_mall_s": time_mall, "time_fd_s": time_fd,
                "fd_selected_bump_frac": fd_out["selected_bump_frac"],
                "fd_status": fd_out["fd_status"],
                "fd_blowup_flag": fd_out["fd_blowup_flag"],
                "fd_blowup_ratio_smallest_to_largest": fd_out["fd_blowup_ratio_smallest_to_largest"],
                "fd_rel_se_selected": fd_out["fd_rel_se_selected"],
                "fd_neighbor_log_gap_selected": fd_out["fd_neighbor_log_gap_selected"],
                "n_steps": config.n_steps, "n_paths_mall": config.n_paths_mall,
                "n_paths_fd": config.n_paths_fd, "n_seeds_fd": config.n_seeds_fd,
            })
            for s in per_seed:
                row = {
                    "H": H, "T": T, "seed": int(s["seed"]), "curvature_fd": s["curvature"], "dk": s["dk"],
                    "price_fit_rmse": s.get("price_fit_rmse", np.nan), "price_fit_cond": s.get("price_fit_cond", np.nan),
                    "iv_fit_rmse": s.get("iv_fit_rmse", np.nan), "iv_fit_cond": s.get("iv_fit_cond", np.nan),
                    "fd_selected_bump_frac": fd_out["selected_bump_frac"], "fd_status": fd_out["fd_status"],
                }
                for key, value in s.items():
                    if key.startswith("iv_") or key.startswith("price_"):
                        row[key] = value
                detail_rows.append(row)

        fit_m = weighted_loglog_fit(config.T_values, mall_vals, mall_ses)
        fit_f = weighted_loglog_fit(config.T_values, fd_vals, fd_ses)
        beta_asian, beta_euro = 2.0 * H + 1.0, 2.0 * H
        pref_th = theoretical_prefactor(H, config.sigma0, config.nu, config.rho)
        pref_m = float(np.exp(fit_m["intercept"])) if np.isfinite(fit_m["intercept"]) else np.nan
        pref_f = float(np.exp(fit_f["intercept"])) if np.isfinite(fit_f["intercept"]) else np.nan
        C_H = get_C_H(H)

        print(f"\n  Malliavin: beta = {fit_m['beta']:.4f} +/- {fit_m['beta_se']:.4f} (R2 = {fit_m['r2']:.6f}), prefactor = {pref_m:.6f}")
        print(f"  FD:        beta = {fit_f['beta']:.4f} +/- {fit_f['beta_se']:.4f} (R2 = {fit_f['r2']:.6f}), prefactor = {pref_f:.6f}")
        print(f"  Theory:    beta = {beta_asian:.2f}, prefactor = {pref_th:.6f}, C(H) = {C_H:.6f}\n")

        exponent_rows.append({
            "H": H, "beta_mall": fit_m["beta"], "se_beta_mall": fit_m["beta_se"], "R2_mall": fit_m["r2"], "intercept_mall": fit_m["intercept"],
            "beta_fd": fit_f["beta"], "se_beta_fd": fit_f["beta_se"], "R2_fd": fit_f["r2"], "intercept_fd": fit_f["intercept"],
            "beta_theory_asian": beta_asian, "beta_theory_euro": beta_euro,
            "beta_error_mall": abs(fit_m["beta"] - beta_asian) if np.isfinite(fit_m["beta"]) else np.nan,
            "beta_error_fd": abs(fit_f["beta"] - beta_asian) if np.isfinite(fit_f["beta"]) else np.nan,
        })
        ratio_m = pref_m / pref_th if np.isfinite(pref_m) and abs(pref_th) > 1e-20 else np.nan
        ratio_f = pref_f / pref_th if np.isfinite(pref_f) and abs(pref_th) > 1e-20 else np.nan
        prefactor_rows.append({"H": H, "C_H": C_H, "pref_mall": pref_m, "pref_fd": pref_f, "pref_theory": pref_th, "pref_ratio_mall": ratio_m, "pref_ratio_fd": ratio_f})
        all_results[H] = {"beta_mall": fit_m["beta"], "se_mall": fit_m["beta_se"], "r2_mall": fit_m["r2"], "beta_fd": fit_f["beta"], "se_fd": fit_f["beta_se"], "r2_fd": fit_f["r2"], "pref_mall": pref_m, "pref_fd": pref_f, "pref_theory": pref_th}

    raw_dir = os.path.join(config.output_dir, "raw")
    summary_dir = os.path.join(config.output_dir, "summaries")
    write_csv(os.path.join(raw_dir, "curvature_raw.csv"), raw_rows)
    write_csv(os.path.join(raw_dir, "scaling_exponents.csv"), exponent_rows)
    write_csv(os.path.join(raw_dir, "prefactors.csv"), prefactor_rows)
    write_csv(os.path.join(raw_dir, "convergence_detail.csv"), detail_rows)
    summary = summarize_validation(config, raw_rows, exponent_rows, prefactor_rows)
    summary["per_H"] = all_results
    write_json(os.path.join(summary_dir, "validation_summary.json"), summary)
    return {"raw_rows": raw_rows, "exponent_rows": exponent_rows, "prefactor_rows": prefactor_rows, "detail_rows": detail_rows, "summary": summary}


def run_european_control(config: ValidationConfig, output_dir: Optional[str] = None) -> Dict[str, Any]:
    outdir = output_dir or config.output_dir
    ensure_output_dirs(outdir)
    cache = DeterministicCache()
    raw_rows: List[Dict[str, Any]] = []
    exponent_rows: List[Dict[str, Any]] = []
    detail_rows: List[Dict[str, Any]] = []

    print("=" * 100)
    print("European Control — Finite-Difference Curvature Validation")
    print("=" * 100)
    for H in config.H_values:
        curv_abs, curv_se = [], []
        print(f"─── H = {H:.2f} ───")
        for T in config.T_values:
            dt = T / config.n_steps
            L = cache.cholesky(H, config.n_steps, dt)
            fd_out = finite_difference_curvature_sweep(
                "european", H, T, config.n_paths_fd, config.n_steps, config.sigma0, config.nu, config.rho,
                n_seeds=config.n_seeds_fd, fd_bump_grid=config.fd_bump_grid, fd_bump_frac_default=config.fd_bump_frac,
                L=L, antithetic=config.antithetic, fd_seed_base=config.fd_seed_base,
                rel_se_max=config.fd_stability_rel_se_max, neighbor_logtol=config.fd_stability_neighbor_logtol,
                strike_shifts=config.fd_smile_shifts,
            )
            curv, se = fd_out["mean_curv"], fd_out["se_curv"]
            curv_abs.append(abs(curv))
            curv_se.append(se)
            raw_rows.append({
                "H": H, "T": T, "curv_fd_euro": curv, "se_fd_euro": se, "abs_curv_fd_euro": abs(curv),
                "log_abs_curv_fd_euro": math.log(abs(curv)) if abs(curv) > 0 else np.nan,
                "beta_theory_euro": 2.0 * H, "fd_selected_bump_frac": fd_out["selected_bump_frac"],
                "fd_status": fd_out["fd_status"], "fd_blowup_flag": fd_out["fd_blowup_flag"],
                "fd_rel_se_selected": fd_out["fd_rel_se_selected"], "fd_neighbor_log_gap_selected": fd_out["fd_neighbor_log_gap_selected"],
                "n_steps": config.n_steps, "n_paths_fd": config.n_paths_fd, "n_seeds_fd": config.n_seeds_fd,
            })
            for s in fd_out["per_seed_selected"]:
                row = {"H": H, "T": T, "seed": int(s["seed"]), "curvature_fd_euro": s["curvature"], "dk": s["dk"]}
                for key, value in s.items():
                    if key.startswith("iv_") or key.startswith("price_") or key.endswith("rmse") or key.endswith("cond"):
                        row[key] = value
                detail_rows.append(row)
            print(f"  T={T:10.5f} | euro FD curv={curv:14.6e} | se={se:10.2e} | theory beta={2*H:.2f}")

        fit = weighted_loglog_fit(config.T_values, curv_abs, curv_se)
        exponent_rows.append({
            "H": H, "beta_fd_euro": fit["beta"], "se_beta_fd_euro": fit["beta_se"], "R2_fd_euro": fit["r2"],
            "intercept_fd_euro": fit["intercept"], "beta_theory_euro": 2.0 * H,
            "beta_error_fd_euro": abs(fit["beta"] - 2.0 * H) if np.isfinite(fit["beta"]) else np.nan,
        })
        print(f"  European FD slope: beta = {fit['beta']:.4f} +/- {fit['beta_se']:.4f} (theory {2*H:.2f}, R2={fit['r2']:.6f})\n")

    raw_dir = os.path.join(outdir, "raw")
    summary_dir = os.path.join(outdir, "summaries")
    diag_dir = os.path.join(outdir, "diagnostics")
    write_csv(os.path.join(raw_dir, "european_curvature_raw.csv"), raw_rows)
    write_csv(os.path.join(raw_dir, "european_scaling_exponents.csv"), exponent_rows)
    write_csv(os.path.join(diag_dir, "european_convergence_detail.csv"), detail_rows)
    summary = summarize_european_control(config, raw_rows, exponent_rows)
    write_json(os.path.join(summary_dir, "european_control_summary.json"), summary)
    return {"raw_rows": raw_rows, "exponent_rows": exponent_rows, "detail_rows": detail_rows, "summary": summary}
