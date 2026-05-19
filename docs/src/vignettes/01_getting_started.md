# Getting Started

End-to-end walk-through of `TVHTE.jl` on a synthetic panel with mixed
cohorts, persistent outcomes, and a policy-reactive covariate.

```@example getting_started
using TVHTE
using Random

Random.seed!(11)
N, T = 1200, 7
cohorts = rand([3.0, 5.0, Inf], N)   # half early, quarter late, quarter never

sim = simulate_tvhte(
    N = N, T = T, t0 = cohorts, J = 3,
    rho_Y = 0.45, rho_delta = 0.55,
    sigma_U = 1.0, sigma_eps = 0.3,
    mu_alpha = 0.0, mu_delta0 = 1.5,
    sigma_alpha = 0.6, sigma_delta0 = 0.5,
    cor_alpha_delta = 0.3,
    beta = 0.4,
    feedback_gamma = [0.1, 0.25, 0.4, 0.3],
    seed = 11
)
size(sim.Y), size(sim.X)
```

## Structural fit

```@example getting_started
fit = tvhte(sim.Y, sim.Y0; t0 = sim.t0, J = sim.J, X = sim.X,
            compute_se = false)
(rho_Y     = round(fit.theta.rho_Y,     digits = 3),
 rho_delta = round(fit.theta.rho_delta, digits = 3),
 beta      = round.(fit.beta, digits = 3),
 mean_delta_path = round.(vec(sum(fit.delta_path, dims = 1) ./ size(fit.delta_path, 1)),
                          digits = 3))
```

## Feedback fit

```@example getting_started
fb = fit_feedback(sim.Y, sim.Y0, sim.X, sim.X0)
(coef = fb.coef, sigma_eta = round(fb.sigma_eta, digits = 3))
```

## Counterfactual decomposition

```@example getting_started
# Treatment shifted earlier (t0* = 2) vs never treated
cf_early = simulate_counterfactual(fit, fb, 2; N_star = 500, seed = 99)
cf_never = simulate_counterfactual(fit, fb, Inf; N_star = 500, seed = 99)

total = vec(sum(cf_early.Y, dims = 1) ./ 500 .-
            sum(cf_never.Y, dims = 1) ./ 500)
round.(total, digits = 3)
```

Total dynamic response (counterfactual `Y[t]` under early treatment
minus never-treated, averaged across simulated units).

The R companion package's `vignettes/illustrative.Rmd` further
decomposes this into direct (through `δ`) vs indirect (through the
covariate feedback) components.
