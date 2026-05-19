# TVHTE.jl

[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://xiangao.github.io/TVHTE.jl/)

Time-varying heterogeneous treatment effects in event studies, in Julia. Companion to the R package [`tvhte`](https://github.com/xiangao/tvhte).

Implements:

- Botosaru, Irene and Laura Liu (2025). "Time-Varying Heterogeneous Treatment Effects in Event Studies." [arXiv:2509.13698](https://arxiv.org/abs/2509.13698).
- Botosaru, Irene and Laura Liu (2026). "Event Studies with Feedback." *AEA Papers and Proceedings* 116: 70–74.

## Why

Standard event-study TWFE regressions assume no residual serial dependence after unit/time fixed effects. When outcomes are persistent, the event-time dummies absorb persistence on top of the causal effect — spurious pre-trends, biased post-treatment estimates.

`TVHTE.jl` fits a dynamic panel with correlated random coefficients on `(α_i, δ_{i0})` and an AR(1) on event-time effects, using two-step semiparametric estimation (QMLE + Gaussian-conjugate empirical Bayes). The (2026) extension factors the likelihood under a homogeneous feedback assumption so direct and indirect dynamic effects can be separately identified.

## Install

```julia
using Pkg
Pkg.add(url = "https://github.com/xiangao/TVHTE.jl")
```

## Documentation & vignettes

Full documentation: **<https://xiangao.github.io/TVHTE.jl/>**

| Page | Description |
|---|---|
| [Home](https://xiangao.github.io/TVHTE.jl/) | Overview, install, motivation |
| [Getting Started](https://xiangao.github.io/TVHTE.jl/dev/vignettes/01_getting_started/) | End-to-end walk-through: simulate → fit → posterior trajectories → feedback → direct/indirect counterfactual decomposition, all with executed output |
| [Reference](https://xiangao.github.io/TVHTE.jl/dev/reference/) | Full API. Each function has its docstring followed by a live `@example` block showing real output (estimates, sizes, etc.) |

## At a glance

```julia
using TVHTE

sim = simulate_tvhte(N = 1500, T = 6, t0 = 3, J = 3,
                     beta = 0.4,
                     feedback_gamma = [0.2, 0.3, 0.5, 0.4], seed = 1)

fit = tvhte(sim.Y, sim.Y0; t0 = sim.t0, J = sim.J, X = sim.X)
fb  = fit_feedback(sim.Y, sim.Y0, sim.X, sim.X0)
cf  = simulate_counterfactual(fit, fb, 5; N_star = 500, seed = 99)
```

The R companion package [`tvhte`](https://github.com/xiangao/tvhte) has a richer vignette including the direct/indirect decomposition plot. The two packages produce matching estimates up to MC noise.
