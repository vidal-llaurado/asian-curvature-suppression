### Asian Curvature Suppression - Draft

Joan Vidal Llauradó, 26 Feb 2026

$$dS_t = \sigma_t\ dW_t, \qquad W = \rho W' + \sqrt{1-\rho^2}W^\perp$$

$$M_t = \frac{1}{T}\int_0^t S_u\ du + \frac{T-t}{T}S_t, \qquad dM_t = \phi_t\ dW_t, \qquad \phi_t = \frac{\sigma_t(T-t)}{T}$$

Fractional Bergomi

$$\sigma_t = \sigma_0\exp\\bigl(\nu\sqrt{2H}\,B_t^H - \tfrac{1}{2}\nu^2 t^{2H}\bigr)$$

$$D_s^{W'}\sigma_u = c_1\\sigma_u(u-s)^{H-1/2}, \qquad c_1 := \tfrac{\nu\sqrt{2H}}{2}$$

$$D_s^{W'}D_r^{W'}\sigma_u = c_1^2\\sigma_u(u-r)^{H-1/2}(u-s)^{H-1/2}$$

Conditional Asian variance

$$v_s^2 = \frac{1}{T-s}\int_s^T \phi_u^2\ du \approx \frac{\sigma_0^2(T-s)^2}{3T^2}$$

---

**First-order weight**

$$\Lambda_s = \phi_s \int_s^T D_s^W\phi_r^2\ dr \sim \frac{2\rho c_1\sigma_0^3}{T^3}(T-s)\int_s^T (T-r)^2(r-s)^{H-1/2}\ dr$$

With $r = s + (T-s)x$:

$$\int_s^T (T-r)^2(r-s)^{H-1/2}\ dr = (T-s)^{H+5/2}\ B(H+\tfrac{1}{2},3)$$

$$\boxed{\Lambda_s \sim \frac{C_1}{T^3}(T-s)^{H+7/2}}$$

---

**Skew verification** ($\beta = H - 1/2$)

Kernel: 

$$\partial_k H \sim v_s^{-3}(T-s)^{-3/2} \sim \frac{T^3}{(T-s)^{9/2}}$$

Integrand: 

$$\partial_k H \cdot \Lambda_s \sim (T-s)^{H-1}$$

$$\int_0^T (T-s)^{H-1}\ ds = \frac{T^H}{H}, \qquad \partial_\sigma B_A \sim T^{1/2}$$

$$\partial_k I_A \sim T^{H-1/2}$$

---

**Second-order weight**

$$\Lambda_s^{(2)} = \phi_s\int_s^T\Bigl[\underbrace{(D_s^W\phi_r)\int_r^T D_r^W\phi_u^2\ du}_{\mathcal{A}} + \underbrace{\phi_r\int_r^T D_s^W D_r^W\phi_u^2\ du}_{\mathcal{B}}\Bigr]dr$$

*Component A.*

$$D_s^W\phi_r = \frac{\rho(T-r)}{T}\ c_1\sigma_r(r-s)^{H-1/2}$$

$$\int_r^T D_r^W\phi_u^2\ du \sim \frac{C\sigma_0^2}{T^2}(T-r)^{H+5/2}$$

$$\mathcal{A}(s) \sim \frac{C_A\sigma_0^4}{T^4}(T-s)\int_s^T (T-r)^{H+7/2}(r-s)^{H-1/2}\ dr$$

The integral: 

$$(T-s)^{2H+4}\ B(H+\tfrac{1}{2},H+\tfrac{9}{2})$$

*Component B.*

$$D_s^W D_r^W\phi_u^2 \sim \frac{4\rho^2 c_1^2\sigma_0^2(T-u)^2}{T^2}(u-r)^{H-1/2}(u-s)^{H-1/2}$$

Inner integral $\sim (T-r)^{H+5/2}(r-s)^{H-1/2}/T^2$, multiplied by $\phi_r \sim (T-r)/T$ and integrated over $r$ yields identical scaling.

$$\boxed{\Lambda_s^{(2)} \sim \frac{C_2}{T^4}(T-s)^{2H+5}}$$

---

**Curvature kernel**

European: $G_E(s) \sim v_s^{-3}(T-s)^{-3/2} \sim (T-s)^{-3/2}$

Asian (double-freezing via $D_s^W v_s^2$): one less singularity

$$\boxed{G_A(s) \sim \frac{1}{v_s^3(T-s)^{1/2}} \sim \frac{T^3}{(T-s)^{7/2}}}$$

---

**Curvature scaling**

$$G_A(s)\cdot\Lambda_s^{(2)} \sim \frac{T^3}{(T-s)^{7/2}} \cdot \frac{(T-s)^{2H+5}}{T^4} = \frac{(T-s)^{2H+3/2}}{T}$$

$$\int_0^T \frac{(T-s)^{2H+3/2}}{T}\ ds = \frac{T^{2H+5/2}}{(2H+\tfrac{5}{2})T} \sim T^{2H+3/2}$$

$$\text{vega}_A \sim T^{1/2} \\Rightarrow\ \boxed{\partial_{kk}^2 I_A \sim T^{2H+1}}$$

---

**European comparison**

| Metric | European | Asian | Net factor |
| --- | --- | --- | --- |
| $\Lambda^{(2)}$ | $(T-s)^{2H+1}$ | $(T-s)^{2H+5}/T^4$ | $(T-s)^4/T^4$ |
| Kernel | $(T-s)^{-3/2}$ | $T^3/(T-s)^{7/2}$ | $T^3/(T-s)^2$ |
| Integrand | $(T-s)^{2H-1/2}$ | $(T-s)^{2H+3/2}/T$ | $(T-s)^2/T$ |
| Integral | $T^{2H+1/2}$ | $T^{2H+3/2}$ | $T^1$ |
| Vega | $T^{1/2}$ | $T^{1/2}$ | $1$ |
| Final | $T^{2H}$ | $T^{2H+1}$ | $T^1$ |

---

**Numerical validation**

| $H$ | $\beta$ (emp) | $2H+1$ | $R^2$ |
| --- | --- | --- | --- |
| 0.1 | 1.187 | 1.2 | 1.0000 |
| 0.2 | 1.394 | 1.4 | 1.0000 |
| 0.3 | 1.600 | 1.6 | 1.0000 |
| 0.4 | 1.802 | 1.8 | 1.0000 |

European-style kernel $G \sim (T-s)^{-3/2}$ yields $\beta \approx 2H$, consistent with $T^{2H}$ scaling. It is still suppressed.

---

$$\text{Skew } \partial_k I_A \sim T^{H-1/2} \quad \text{(same singularity as European)}$$

$$\text{Curvature } \partial_{kk}^2 I_A \sim T^{2H+1} \quad \text{(European: } T^{2H}\text{; suppressed by } T^1\text{)}$$

The Asian diffusion coefficient $\phi_t$ contributes four powers of $(T-s)/T$ at second order. The kernel absorbs $T^3/(T-s)^2$ of this, leaving a net regularizing integrand factor of $(T-s)^2/T$. Upon integration, this yields a strict $T^1$ suppression relative to the European case.
