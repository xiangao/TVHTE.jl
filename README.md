# TVHTE.jl

[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://xiangao.github.io/TVHTE.jl/)

`TVHTE.jl` estimates time-varying heterogeneous treatment effects in event
studies. It is the Julia companion to the R package
[`tvhte`](https://github.com/xiangao/tvhte).

Implements:

- Botosaru, Irene and Laura Liu (2025). "Time-Varying Heterogeneous Treatment Effects in Event Studies." [arXiv:2509.13698](https://arxiv.org/abs/2509.13698).
- Botosaru, Irene and Laura Liu (2026). "Event Studies with Feedback." *AEA Papers and Proceedings* 116: 70–74.

## Why

Standard event-study TWFE regressions leave little room for residual serial
dependence. When outcomes are persistent, event-time dummies can pick up
persistence as well as the treatment effect. The result can look like
pre-trends or biased post-treatment effects.

`TVHTE.jl` fits the dynamic panel model with correlated random coefficients and
an AR(1) structure on event-time effects. The feedback extension uses the
homogeneous-feedback restriction from Botosaru and Liu (2026) to separate
direct and indirect dynamic effects.

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

The R companion package [`tvhte`](https://github.com/xiangao/tvhte) has the
longer vignette with the direct/indirect decomposition plot. The two packages
are tested on the same simulated designs.
