"""
Monte Carlo: Asian curvature suppression

Two independent methods:
  - Pathwise second-order Malliavin weight integration
  - Finite-difference implied volatility curvature
"""

import numpy as np
from scipy.stats import norm, linregress
from scipy.optimize import brentq
import time
import warnings
warnings.filterwarnings('ignore')



# fBm generation via Cholesky

def generate_fBm_cholesky(H, n_steps, dt, n_paths):
    def fbm_cov(s, t, H):
        return 0.5 * (abs(s)**(2*H) + abs(t)**(2*H) - abs(t - s)**(2*H))

    incr_cov = np.zeros((n_steps, n_steps))
    for i in range(n_steps):
        for j in range(n_steps):
            ti0, ti1 = i * dt, (i + 1) * dt
            tj0, tj1 = j * dt, (j + 1) * dt
            incr_cov[i, j] = (fbm_cov(ti1, tj1, H) - fbm_cov(ti1, tj0, H)
                              - fbm_cov(ti0, tj1, H) + fbm_cov(ti0, tj0, H))

    incr_cov += 1e-14 * np.eye(n_steps)
    L = np.linalg.cholesky(incr_cov)
    Z = np.random.randn(n_paths, n_steps)
    return Z @ L.T


def simulate_fbergomi(H, T, n_paths, n_steps, sigma0, nu, rho, antithetic=True):
    dt = T / n_steps
    t_grid = np.linspace(0, T, n_steps + 1)

    half = n_paths // 2 if antithetic else n_paths

    # fBm increments for W'
    dBH = generate_fBm_cholesky(H, n_steps, dt, half)

    if antithetic:
        dBH = np.vstack([dBH, -dBH])
        n_eff = 2 * half
    else:
        n_eff = half

    # Construct fBm path
    BH = np.zeros((n_eff, n_steps + 1))
    BH[:, 1:] = np.cumsum(dBH, axis=1)

    # Volatility process (log-normal)
    V = sigma0**2 * np.exp(
        nu * np.sqrt(2 * H) * BH - 0.5 * nu**2 * t_grid**(2*H)
    )
    sigma = np.sqrt(np.maximum(V, 1e-20))

    # Asian diffusion coefficient
    phi = np.zeros_like(sigma)
    for i in range(n_steps + 1):
        phi[:, i] = sigma[:, i] * (T - t_grid[i]) / T

    return sigma, phi, t_grid, dBH



# Pathwise Malliavin weight method

def compute_Lambda2_pathwise(sigma, phi, T, n_steps, dt, H, nu, rho):
    n_paths = sigma.shape[0]
    c1 = nu * np.sqrt(2 * H) / 2.0

    Lambda2 = np.zeros((n_paths, n_steps))

    for s_idx in range(n_steps - 2):
        s = (s_idx + 0.5) * dt  # midpoint rule
        phi_s = phi[:, s_idx]
        integral_r = np.zeros(n_paths)

        for r_idx in range(s_idx + 1, n_steps - 1):
            r = (r_idx + 0.5) * dt
            phi_r = phi[:, r_idx]
            sigma_r = sigma[:, r_idx]

            # D_s^W phi_r = rho * (T-r)/T * sigma_r * c1 * (r-s)^{H-1/2}
            rs_kernel = max(r - s, dt * 0.1) ** (H - 0.5)
            D_s_W_phi_r = rho * ((T - r) / T) * sigma_r * c1 * rs_kernel

            # Inner u-integral
            u_indices = np.arange(r_idx + 1, n_steps)
            if len(u_indices) == 0:
                continue

            u_vals = (u_indices + 0.5) * dt
            sigma_u = sigma[:, u_indices]  # (n_paths, n_u)
            weight_u = (T - u_vals) ** 2 / T ** 2  # (n_u,)

            ur_kernel = np.maximum(u_vals - r, dt * 0.1) ** (H - 0.5)
            us_kernel = np.maximum(u_vals - s, dt * 0.1) ** (H - 0.5)

            # Component A inner
            integrand_A = 2 * rho * weight_u * (sigma_u ** 2) * c1 * ur_kernel
            inner_A = np.sum(integrand_A, axis=1) * dt

            # Component B inner
            integrand_B = 4 * rho**2 * weight_u * (sigma_u**2) * c1**2 * ur_kernel * us_kernel
            inner_B = np.sum(integrand_B, axis=1) * dt

            comp_A = D_s_W_phi_r * inner_A
            comp_B = phi_r * inner_B

            integral_r += (comp_A + comp_B) * dt

        Lambda2[:, s_idx] = phi_s * integral_r

    return Lambda2


def compute_curvature_correction(Lambda2, sigma, phi, T, n_steps, dt, sigma0):
    n_paths = Lambda2.shape[0]
    integrand = np.zeros((n_paths, n_steps))

    for s_idx in range(n_steps - 2):
        s = (s_idx + 0.5) * dt
        tau = T - s
        if tau < dt * 0.5:
            continue

        # Conditional vol
        v_s = sigma0 * tau / (T * np.sqrt(3))

        # Asian curvature kernel
        G_s = 1.0 / (v_s**3 * tau**0.5)

        integrand[:, s_idx] = G_s * Lambda2[:, s_idx]

    integral = np.sum(integrand, axis=1) * dt

    # Asian Bachelier Vega
    vega_A = np.sqrt(T / (2 * np.pi)) * sigma0 / np.sqrt(3)

    curvature_samples = integral / vega_A
    mean_corr = np.mean(curvature_samples)
    se_corr = np.std(curvature_samples) / np.sqrt(n_paths)

    return mean_corr, se_corr, curvature_samples



# Finite-difference implied volatility curvature

def bachelier_price(x, k, sigma, T):
    if sigma * np.sqrt(T) < 1e-15:
        return max(x - k, 0.0)
    d = (x - k) / (sigma * np.sqrt(T))
    return sigma * np.sqrt(T) * (d * norm.cdf(d) + norm.pdf(d))


def bachelier_iv_inversion(price, x, k, T, iv_guess):
    if price < 1e-15:
        return iv_guess
    try:
        def obj(sig):
            return bachelier_price(x, k, sig, T) - price
        lo, hi = 1e-6, 10 * iv_guess
        while obj(hi) < 0 and hi < 1e6:
            hi *= 2
        return brentq(obj, lo, hi, xtol=1e-12)
    except:
        return iv_guess


def finite_difference_curvature(H, T, n_paths, n_steps, sigma0, nu, rho, n_seeds=5):
    dt = T / n_steps
    v_asian_0 = sigma0 / np.sqrt(3)  # ATM Asian IV at leading order
    dk = v_asian_0 * np.sqrt(T) * 0.10  # 10% of ATM Asian stdev

    curvatures = []

    for seed in range(n_seeds):
        np.random.seed(1000 * seed + 42)

        # Generate volatility paths
        sigma, phi, t_grid, _ = simulate_fbergomi(
            H, T, n_paths, n_steps, sigma0, nu, rho, antithetic=True
        )

        # Generate asset noise
        dW_prime = np.random.randn(sigma.shape[0], n_steps) * np.sqrt(dt)
        dW_perp = np.random.randn(sigma.shape[0], n_steps) * np.sqrt(dt)
        dW = rho * dW_prime + np.sqrt(1 - rho**2) * dW_perp

        # Construct asset paths
        S = np.zeros((sigma.shape[0], n_steps + 1))
        for i in range(n_steps):
            S[:, i + 1] = S[:, i] + sigma[:, i] * dW[:, i]

        # Asian average (trapezoidal)
        weights_trap = np.ones(n_steps + 1)
        weights_trap[0] = weights_trap[-1] = 0.5
        A_T = np.sum(S * weights_trap, axis=1) * dt / T

        # Price at three strikes
        ivs = {}
        for shift in [-1, 0, 1]:
            k = shift * dk
            payoff = np.maximum(A_T - k, 0)
            price = np.mean(payoff)
            iv = bachelier_iv_inversion(price, 0.0, k, T, v_asian_0)
            ivs[shift] = iv

        curv = (ivs[1] - 2 * ivs[0] + ivs[-1]) / dk**2
        curvatures.append(curv)

    mean_curv = np.mean(curvatures)
    se_curv = np.std(curvatures) / np.sqrt(len(curvatures))
    return mean_curv, se_curv



# Main validation

def run_full_validation():
    H_values = [0.1, 0.2, 0.3, 0.4]
    T_values = [0.1, 0.05, 0.025, 0.0125, 0.00625]

    # Simulation parameters
    n_paths = 8000
    n_steps = 64
    sigma0, nu, rho = 0.3, 0.5, -0.7

    print(f"Paths: {n_paths} (antithetic); Steps: {n_steps}")
    print(f"Sigma0 = {sigma0}, Nu = {nu}, Rho = {rho}")

    all_results = {}

    for H in H_values:
        print(f"  H = {H}")
        print(f"{'T':>10} | {'|Malliavin|':>14} | {'SE_M':>10} | "
              f"{'|FD curv|':>14} | {'SE_FD':>10} | {'logT':>8}")
        print("-" * 80)

        mall_corrections = []
        fd_corrections = []
        log_Ts = []

        for T in T_values:
            dt = T / n_steps
            np.random.seed(42)

            # Malliavin pathwise
            sigma, phi, t_grid, _ = simulate_fbergomi(
                H, T, n_paths, n_steps, sigma0, nu, rho
            )
            Lambda2 = compute_Lambda2_pathwise(sigma, phi, T, n_steps, dt, H, nu, rho)
            mall_mean, mall_se, _ = compute_curvature_correction(
                Lambda2, sigma, phi, T, n_steps, dt, sigma0
            )
            abs_mall = abs(mall_mean)

            # Finite difference
            fd_mean, fd_se = finite_difference_curvature(
                H, T, n_paths // 2, n_steps, sigma0, nu, rho, n_seeds=3
            )
            abs_fd = abs(fd_mean)

            mall_corrections.append(abs_mall)
            fd_corrections.append(abs_fd)
            log_Ts.append(np.log(T))

            print(f"{T:10.4f} | {abs_mall:14.6e} | {mall_se:10.2e} | "
                  f"{abs_fd:14.6e} | {fd_se:10.2e} | {np.log(T):8.2f}")

        # Fit scaling exponents
        log_mall = np.log(np.array(mall_corrections))
        log_fd = np.log(np.maximum(np.array(fd_corrections), 1e-20))
        log_T = np.array(log_Ts)

        # Malliavin method fit
        valid_m = np.isfinite(log_mall)
        if np.sum(valid_m) >= 3:
            slope_m, intercept_m, r_m, p_m, se_m = linregress(log_T[valid_m], log_mall[valid_m])
        else:
            slope_m, r_m, se_m = np.nan, np.nan, np.nan

        # FD method fit
        valid_f = np.isfinite(log_fd) & (np.array(fd_corrections) > 1e-15)
        if np.sum(valid_f) >= 3:
            slope_f, intercept_f, r_f, p_f, se_f = linregress(log_T[valid_f], log_fd[valid_f])
        else:
            slope_f, r_f, se_f = np.nan, np.nan, np.nan

        theoretical_asian = 2 * H + 1
        theoretical_euro = 2 * H - 1

        print(f"\n  Malliavin method (Asian kernel τ^{{-1/2}}):")
        print(f"    beta_emp = {slope_m:.3f} ± {se_m:.3f}  (R2 = {r_m**2:.4f})")
        print(f"  Finite difference method:")
        print(f"    beta_emp = {slope_f:.3f} ± {se_f:.3f}  (R2 = {r_f**2:.4f})")
        print(f"  Theoretical:")
        print(f"    Asian (predicted): beta = 2H + 1 = {theoretical_asian:.1f}")
        print(f"    European (known):  beta = 2H = {theoretical_euro:.1f}")

        if slope_m > 0:
            match = abs(slope_m - theoretical_asian) < 0.15
            if match:
                print(f"  Matches 2H + 1: beta = {slope_m:.3f} ≈ {theoretical_asian:.1f}")
            else:
                print(f"  Supressed (beta = {slope_m:.3f} > 0) but offset from 2H + 1 = {theoretical_asian:.1f}")
        else:
            print(f"  beta_M = {slope_m:.2f} ≤ 0")

        all_results[H] = {
            'beta_mall': slope_m, 'se_mall': se_m, 'r2_mall': r_m**2,
            'beta_fd': slope_f, 'se_fd': se_f, 'r2_fd': r_f**2 if not np.isnan(r_f) else np.nan,
            'theoretical': theoretical_asian,
            'mall_corrections': mall_corrections,
            'fd_corrections': fd_corrections,
        }


    # Summary table
    print(f"{'H':>5} | {'beta_Mall':>10} | {'beta_FD':>10} | {'beta_theory':>10} | "
          f"{'beta_Euro':>10} | {'R2_M':>8} | {'Suppressed':>12}")
    print("-" * 82)
    for H in H_values:
        r = all_results[H]
        supp = "Yes" if r['beta_mall'] > 0 else "No"
        print(f"{H:5.1f} | {r['beta_mall']:10.3f} | {r['beta_fd']:10.3f} | "
              f"{r['theoretical']:10.1f} | {2*H-1:10.1f} | "
              f"{r['r2_mall']:8.4f} | {supp:>12}")

    for H in H_values:
        r = all_results[H]
        t_stat = r['beta_mall'] / r['se_mall'] if r['se_mall'] > 0 else np.inf
        print(f"  H={H}: beta = {r['beta_mall']:.3f} ± {r['se_mall']:.3f}, "
              f"t-stat = {t_stat:.1f}, "
              f"{'Reject beta ≤ 0' if t_stat > 2 else 'Cannot reject beta ≤ 0'} at 95%")

    return all_results


if __name__ == "__main__":
    results = run_full_validation()
