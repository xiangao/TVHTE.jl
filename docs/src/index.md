# TVHTE.jl

Time-varying heterogeneous treatment effects in event studies. Julia
port of the R companion [`tvhte`](https://github.com/xiangao/tvhte).

Implements:

- **Botosaru, Irene and Laura Liu (2025).** "Time-Varying Heterogeneous
  Treatment Effects in Event Studies."
  [arXiv:2509.13698](https://arxiv.org/abs/2509.13698).
- **Botosaru, Irene and Laura Liu (2026).** "Event Studies with
  Feedback." *AEA Papers and Proceedings* 116: 70–74.

## Why

Standard event-study TWFE regressions implicitly assume the residual
has no serial dependence after unit/time fixed effects. When outcomes
are persistent — earnings, employment, consumption, anything with
habits or adjustment costs — the event-time dummies absorb both the
true causal path **and** residual dynamics, producing spurious
pre-trends and biased post-treatment estimates.

`TVHTE.jl` fits a dynamic panel with **correlated random coefficients
on `(α_i, δ_{i0})`** and an **AR(1) on the event-time effects**:

```
Y_it    = ρ_Y Y_{i,t-1} + α_i + Σ_j D_{it}^j δ_{ij} + X_{it}'β + U_it
δ_{ij}  = ρ_δ δ_{i,j-1} + ε_{ij},   j ≥ 1
```

A two-step semiparametric estimator: QMLE for common parameters
(integrates `λ_i = (α_i, δ_{i0})` analytically under a Gaussian
working assumption); Gaussian-conjugate empirical Bayes for unit-level
posterior trajectories `{δ_{i,j}}_{j=0}^J`.

The Botosaru-Liu (2026) extension models the covariate `X` as
endogenous to past outcomes via a homogeneous feedback process; the
likelihood factors into a structural piece and a feedback piece that
are separately identified.

## Three estimators

- [`tvhte`](@ref) — main two-step estimator. Supports staggered
  adoption (`t0` as a vector with `Inf` for never-treated),
  strictly-exogenous covariates `X`, and optional standard errors.
- [`fit_feedback`](@ref) — pooled OLS for the AR(1) feedback process.
- [`simulate_counterfactual`](@ref) — Botosaru-Liu (2026) Algorithm 1:
  joint counterfactual `(Y*, X*)` under alternative treatment timing,
  decomposing dynamic responses into direct vs indirect components.

## Install

```julia
using Pkg
Pkg.add(url = "https://github.com/xiangao/TVHTE.jl")
```

## At a glance

```julia
using TVHTE

sim = simulate_tvhte(N = 1500, T = 6, t0 = 3, J = 3,
                     rho_Y = 0.5, rho_delta = 0.7,
                     beta = 0.4,
                     feedback_gamma = [0.2, 0.3, 0.5, 0.4], seed = 1)

fit = tvhte(sim.Y, sim.Y0; t0 = sim.t0, J = sim.J, X = sim.X)
fb  = fit_feedback(sim.Y, sim.Y0, sim.X, sim.X0)

cf  = simulate_counterfactual(fit, fb, 5; N_star = 500, seed = 99)
```

See the **Vignettes** in the sidebar for end-to-end examples and the
[Reference](reference.md) page for the API.

The R companion [`tvhte`](https://github.com/xiangao/tvhte) follows
the same API and is tested against the same DGPs; estimates match to
MC noise.
