# Getting Started

End-to-end walk-through of `TVHTE.jl` on a synthetic panel with mixed
cohorts (including never-treated units), persistent outcomes,
unit-level heterogeneity in dynamic responses, and a policy-reactive
covariate that adjusts in response to past outcomes. This mirrors the
R companion's vignette, with the same DGP and the same numeric
decomposition.

The DGP here is synthetic only because the Botosaru-Liu (2025)
county-unemployment example doesn't ship replication data. On a real
panel you swap the simulated `Y`, `Y0`, `X`, `X0`, `t0` for your own
arrays of the same shape.

## Setup and DGP

```@example getting_started
using TVHTE
using Random
using Statistics

Random.seed!(11)
N, T = 1500, 7
# Mixed cohort: half treated at t=3, quarter at t=5, quarter never.
cohorts = rand([3.0, 5.0, Inf], N)
sum(cohorts .== Inf), sum(cohorts .< Inf)
```

```@example getting_started
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
size(sim.Y), size(sim.X), length(sim.t0)
```

The covariate `X_{it}` evolves each period in response to lagged `Y`
and lagged `X` — the homogeneous-feedback DGP of Botosaru-Liu (2026).

## Structural fit

```@example getting_started
fit = tvhte(sim.Y, sim.Y0; t0 = sim.t0, J = sim.J, X = sim.X,
            compute_se = false)

(rho_Y     = round(fit.theta.rho_Y,     digits = 3),
 rho_delta = round(fit.theta.rho_delta, digits = 3),
 sigma_U   = round(sqrt(fit.theta.sigma_U2),   digits = 3),
 sigma_eps = round(sqrt(fit.theta.sigma_eps2), digits = 3),
 beta      = round.(fit.beta, digits = 3))
```

True values were `rho_Y = 0.45`, `rho_delta = 0.55`, `beta = 0.4`.

The prior moments of the correlated random coefficients
`(alpha_i, delta_{i0})`:

```@example getting_started
(mu = (round(fit.prior.mu_alpha, digits = 3),
       round(fit.prior.mu_delta0, digits = 3)),
 sd = (round(sqrt(fit.prior.sigma_alpha2),   digits = 3),
       round(sqrt(fit.prior.sigma_delta0_2), digits = 3)),
 cor = round(fit.prior.cov_alpha_delta /
             sqrt(fit.prior.sigma_alpha2 * fit.prior.sigma_delta0_2),
             digits = 3))
```

True: `mu = (0, 1.5)`, `sd = (0.6, 0.5)`, `cor = 0.3`.

## Posterior event-time trajectories

The empirical-Bayes step gives a per-unit `(alpha_hat, delta0_hat)`
plus the full event-time path `delta_{i,j}` for `j = 0:J`. Averaged
across treated units, this is the kind of event-study path an applied
researcher would plot:

```@example getting_started
treated = findall(isfinite, sim.t0)
mean_path = vec(mean(view(fit.delta_path, treated, :), dims = 1))
round.(mean_path, digits = 3)
```

The path decays at roughly the estimated `rho_delta`.

The unit-level posteriors fan out around this average because
`(alpha_i, delta_{i0})` varies across units — this is what
distinguishes TV-HTE from the single-`delta_j` standard event study.

```@example getting_started
sub = treated[1:5]
round.(fit.delta_path[sub, :], digits = 3)
```

## Feedback fit

```@example getting_started
fb = fit_feedback(sim.Y, sim.Y0, sim.X, sim.X0)
(coef = (intercept = round(fb.coef.intercept, digits = 3),
         gamma_Y   = round(fb.coef.gamma_Y,   digits = 3),
         gamma_X   = round(fb.coef.gamma_X,   digits = 3)),
 sigma_eta = round(fb.sigma_eta, digits = 3))
```

The covariate is modelled as an AR(1) with one lag of `Y`. True
coefficients were `(0.1, 0.25, 0.4)` with `sigma_eta = 0.3`.

Under the BL 2026 factorisation, the structural piece (fit by
`tvhte()`) and the feedback piece (fit by `fit_feedback()`) are
estimated separately — neither needs a parametric specification of the
other.

## Counterfactual decomposition

`simulate_counterfactual` implements BL 2026 Algorithm 1: it draws
latent types from the estimated prior, evolves event-time effects via
the AR(1), and propagates `(Y, X)` jointly using the estimated
feedback process. Output is a joint counterfactual path
`(Y_i^{T,*}, X_i^{T,*})`, separating the **direct** effect (through
`delta_{i,j}`) from the **indirect** effect (through `X_t * beta`).

```@example getting_started
# Counterfactual 1: same units, treatment shifted earlier (t0 = 2)
cf_early = simulate_counterfactual(fit, fb, 2; N_star = 500, seed = 99)
# Counterfactual 2: never treated
cf_never = simulate_counterfactual(fit, fb, Inf; N_star = 500, seed = 99)

total_effect = vec(mean(cf_early.Y, dims = 1)) .-
               vec(mean(cf_never.Y, dims = 1))
round.(total_effect, digits = 3)
```

The total dynamic response: counterfactual `Y[t]` averaged across
simulated units, under early treatment minus under no treatment.

To isolate the **direct** component, run the same counterfactual with
the feedback coefficients zeroed out (covariate held to its noise
path, no response to lagged outcomes):

```@example getting_started
fb_frozen = (coef = (intercept = fb.coef.intercept,
                     gamma_Y   = 0.0,
                     gamma_X   = 0.0),
             sigma_eta = fb.sigma_eta)

cf_early_direct = simulate_counterfactual(fit, fb_frozen, 2;
                                           N_star = 500, seed = 99)
cf_never_direct = simulate_counterfactual(fit, fb_frozen, Inf;
                                           N_star = 500, seed = 99)
direct_effect = vec(mean(cf_early_direct.Y, dims = 1)) .-
                vec(mean(cf_never_direct.Y, dims = 1))
indirect_effect = total_effect .- direct_effect

[round.(total_effect,    digits = 3)';
 round.(direct_effect,   digits = 3)';
 round.(indirect_effect, digits = 3)']
```

Rows: total response, direct (through `delta`), indirect (through
covariate feedback). In a minimum-wage application, the direct row is
the wage response holding firm input demand fixed; the indirect row is
the additional response that flows through firms re-optimising hours
and employment composition.

## On real data

The workflow is the same:

1. Arrange your panel as an `N × T` outcome matrix `Y`, baseline
   vector `Y0`, optional `N × T × K` covariate array `X`, baseline
   `X0`, and length-`N` treatment-cohort vector `t0` (use `Inf` for
   never-treated units).
2. `tvhte(Y, Y0; t0, J, X)` for the structural fit.
3. `fit_feedback(Y, Y0, X, X0)` for the feedback fit.
4. `simulate_counterfactual(fit, fit_fb, t0_star)` for the
   counterfactual decomposition.

The QMLE step under the Gaussian working prior is consistent for the
common parameters even when the true prior on
`(alpha_i, delta_{i0})` is non-Gaussian (BL 2025). The Gaussian-
conjugate empirical-Bayes step achieves ratio optimality. The
homogeneous-feedback factorisation delivers separately-identified
direct and indirect effects (BL 2026).

## References

- Botosaru, Irene and Laura Liu (2025). "Time-Varying Heterogeneous
  Treatment Effects in Event Studies."
  [arXiv:2509.13698](https://arxiv.org/abs/2509.13698).
- Botosaru, Irene and Laura Liu (2026). "Event Studies with
  Feedback." *AEA Papers and Proceedings* 116: 70–74.
