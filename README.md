## Asian Curvature Suppression under Rough Volatility

Numerical validation of the short-maturity scaling of ATM implied-volatility curvature for Asian options under rough Bachelier dynamics, using Malliavin methods and controlled Monte Carlo experiments.

---

#### Target result

The ATM Bachelier implied-volatility curvature of an Asian option satisfies

$$\partial_{kk}^2 I_A(T, 0) \sim C(H, \sigma_0, \nu, \rho)\, T^{2H+1}, \qquad T \downarrow 0,$$

which is one power of $T$ slower than the European curvature $\partial_{kk}^2 I_E(T, 0) \sim \tilde{C}(H, \sigma_0, \nu, \rho)\, T^{2H}$.

This suppression is a consequence of time-averaging. The Asian payoff applies an order-one integration operator to the variance path, raising its Sobolev regularity by $+1$ and shifting the curvature exponent accordingly.

#### Model

Bachelier–rough Bergomi dynamics under $\mathbb{Q}$:

$$dS_t = \sigma_t\, dW_t, \qquad \sigma_t^2 = \sigma_0^2 \exp\!\Big(\nu \sqrt{2H}\, B_t^H - \tfrac{1}{2}\nu^2 t^{2H}\Big),$$

where $B^H$ is Riemann–Liouville fractional Brownian motion with Hurst parameter $H \in (0, \tfrac{1}{2})$, and $W$ is correlated with $B^H$ via $\rho \in (-1, 0)$.

The Asian payoff depends on the time-average $A_T = T^{-1}\!\int_0^T S_t\, dt$. The quantity of interest is the ATM curvature $\kappa_A(T) := \partial_{kk}^2 I_A(T, k)\big|_{k=0}$.

#### Origin of the extra power of $T$

The effective asymptotic Asian diffusion coefficient is $\varphi_t = \sigma_t (T - t) / T$. The curvature involves quadratic functionals of $\varphi$: relative to the European case, each occurrence of $\sigma_t$ is replaced by $\sigma_t(T - t)/T$, introducing additional time weights.

At second order, four powers of $(T - t)/T$ appear. Two are absorbed by the time integration over $[0, T]$, and two by the structure of the curvature kernel. The net balance is a factor of $T^1$.

This is the semigroup composition $\mathcal{V}_1 \circ \mathcal{V}_{H+1/2} = c\, \mathcal{V}_{H+3/2}$. The averaging operator raises the effective regularity by exactly $+1$, independently of $H$, producing the shift from $2H$ to $2H + 1$ in the curvature exponent.

#### Prefactor

The leading-order prefactor is

$$C(H, \sigma_0, \nu, \rho) = \frac{3\sqrt{6\pi}\, \rho^2 \nu^2 H\, \sigma_0}{2H + \tfrac{5}{2}} \cdot \mathcal{C}(H),$$

where the geometric constant is

$$\mathcal{C}(H) = B\!\left(H + \tfrac{1}{2},\, 3\right) B\!\left(H + \tfrac{1}{2},\, H + \tfrac{9}{2}\right) + 2K(H).$$

The term $K(H)$ is a two-dimensional rough-kernel contribution that does not factorise into Beta functions for general $H$:

$$K(H) = \int_0^1 \!\int_0^1 (1 - \xi)^{H+7/2}\, t^{H-1/2}\, (1-t)^2\, \bigl(\xi + (1-\xi)t\bigr)^{H-1/2}\, dt\, d\xi.$$

The non-factorisation arises because the two fractional powers $(u - r)^{H-1/2}(u - s)^{H-1/2}$ in the second Malliavin derivative create simultaneous dependence on both reference points.

At $H = \tfrac{1}{2}$: $K(\tfrac{1}{2}) = \tfrac{1}{15}$ and $\mathcal{C}(\tfrac{1}{2}) = \tfrac{1}{5}$.

In the implementation, $\mathcal{C}(H)$ is evaluated via high-order quadrature, the full prefactor is computed explicitly, and Monte Carlo estimates of $\kappa_A(T)$ are compared against $C(H, \sigma_0, \nu, \rho)\, T^{2H+1}$.

### Numerical methodology

#### Malliavin estimator

The curvature is computed via a pathwise second-order Malliavin weight. This evaluates the asymptotic contribution directly, avoids numerical differentiation, and remains stable across the maturity range tested.

This is the primary validation method.

#### Finite-difference implied-volatility estimator

A secondary estimator based on local strike grids, price-smile fitting, implied-vol inversion, and quadratic curvature extraction. Even after stabilisation (bump sweeps, local fits), this approach becomes unreliable as $T \to 0$.

The issue is that the curvature signal scales as $T^{2H+1} \to 0$, Monte Carlo noise does not decay at the same rate, implied-vol inversion amplifies the noise, and second derivatives are ill-conditioned. The FD route is therefore included as a diagnostic benchmark, not as primary evidence.

#### European control

A parallel experiment estimates the European curvature exponent $2H$ using the same finite-difference pipeline. This provides a baseline and a direct numerical comparison with the Asian exponent $2H + 1$.

### Repository structure


- `asian_validation_common.py` Numerical engine and shared utilities
- `run_validation.py` Main Asian validation run
- `run_european_control.py` European control experiment
- `run_robustness.py` Step-count and path-budget robustness sweeps
- `make_plots.py`
- `run_all.sh`

Outputs are written under `outputs/` in subdirectories `raw/`, `summaries/`, `figures/`, and `diagnostics/`.

### Running the validation

#### Smoke test

Verifies the pipeline, file generation, and plotting:

```bash
python run_validation.py --mode smoke --output-dir outputs
python run_european_control.py --mode smoke --output-dir outputs
python make_plots.py
```

Smoke mode uses `H ∈ {0.10, 0.30}`, four maturities, 1000 Malliavin paths, and 128 time steps.

#### Complete validation

```bash
python run_validation.py \
  --mode full \
  --output-dir outputs \
  --n-paths-mall 8000 \
  --n-paths-fd 4000 \
  --n-steps 256 \
  --n-seeds-fd 7 \
  --fd-bump-grid 0.1,0.2,0.4

python run_european_control.py \
  --mode full \
  --output-dir outputs \
  --n-paths-fd 6000 \
  --n-steps 256 \
  --n-seeds-fd 7 \
  --fd-bump-grid 0.1,0.2,0.4

python make_plots.py
```


#### Baseline parameters

| Parameter | Value |
|---|---|
| $H$ | 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40 |
| $T$ | 0.1, 0.05, 0.025, 0.0125, 0.00625 |
| $\sigma_0$ | 0.3 |
| $\nu$ | 0.5 |
| $\rho$ | −0.7 |

### Outputs

**Raw data.** `curvature_raw.csv`, `scaling_exponents.csv`, `prefactors.csv`, `european_curvature_raw.csv`, `european_scaling_exponents.csv`.

**Summaries.** `validation_summary.json`, `european_control_summary.json`.

**Figures.** log–log scaling plots, normalised curvature, exponent errors, Malliavin and FD comparison, European and Asian comparison.

### Interpreting the results

The primary objects to inspect are:

- The fitted Malliavin slope versus the target $2H + 1$.
- The normalised curvature plateau $\kappa_{\mathrm{mall}} / T^{2H+1}$, which should stabilise near the theoretical prefactor.
- The European FD slope versus $2H$, confirming the one-power gap.
- The cross-method discrepancy plot, which documents the FD estimator's degradation at short maturities.

If the FD method fails to recover the Asian exponent while the Malliavin estimator remains stable, that is consistent with the known ill-conditioning of extracting implied-volatility curvature from Monte Carlo prices in this regime.
