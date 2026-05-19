# API Reference

```@meta
CurrentModule = TVHTE
```

Each function's docstring describes the API; an `@example` block
immediately after shows it running with real output. For an end-to-end
walk-through, see the
[Getting Started](vignettes/01_getting_started.md) vignette.

## Main estimator

```@docs
tvhte
```

### Example

```@example tvhte
using TVHTE
using Random
Random.seed!(1)

sim = simulate_tvhte(N = 300, T = 6, t0 = 3, J = 3,
                     rho_Y = 0.5, rho_delta = 0.7, seed = 1)
fit = tvhte(sim.Y, sim.Y0; t0 = sim.t0, J = sim.J, compute_se = false)
(rho_Y = round(fit.theta.rho_Y, digits = 3),
 rho_delta = round(fit.theta.rho_delta, digits = 3))
```

## Covariate feedback (Botosaru-Liu 2026)

```@docs
fit_feedback
```

### Example

```@example fb
using TVHTE
using Random
Random.seed!(1)

sim = simulate_tvhte(N = 300, T = 6, t0 = 3, J = 2,
                     beta = 0.4,
                     feedback_gamma = [0.2, 0.3, 0.5, 0.4], seed = 1)
fb = fit_feedback(sim.Y, sim.Y0, sim.X, sim.X0)
(coef = (intercept = round(fb.coef.intercept, digits = 3),
         gamma_Y   = round(fb.coef.gamma_Y,   digits = 3),
         gamma_X   = round(fb.coef.gamma_X,   digits = 3)),
 sigma_eta = round(fb.sigma_eta, digits = 3))
```

```@docs
simulate_counterfactual
```

### Example

```@example cf
using TVHTE
using Random
using Statistics
Random.seed!(1)

sim = simulate_tvhte(N = 300, T = 6, t0 = 3, J = 2, beta = 0.4,
                     feedback_gamma = [0.2, 0.3, 0.5, 0.4], seed = 1)
fit = tvhte(sim.Y, sim.Y0; t0 = sim.t0, J = sim.J, X = sim.X,
            compute_se = false)
fb  = fit_feedback(sim.Y, sim.Y0, sim.X, sim.X0)

# Counterfactual: shift treatment timing from 3 to 5
cf = simulate_counterfactual(fit, fb, 5; N_star = 200, seed = 99)
(size_Y = size(cf.Y),
 mean_Y_path = round.(vec(mean(cf.Y, dims = 1)), digits = 3))
```

## Simulation

```@docs
simulate_tvhte
```

### Example

```@example sim
using TVHTE

# Common adoption
sim = simulate_tvhte(N = 200, T = 5, t0 = 3, J = 2,
                     rho_Y = 0.4, rho_delta = 0.6, seed = 1)
(size_Y = size(sim.Y), first_unit_Y = round.(sim.Y[1, :], digits = 3))
```
