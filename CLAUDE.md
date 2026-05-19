# TVHTE.jl — project notes for Claude

Julia port of [`tvhte`](https://github.com/xiangao/tvhte). Botosaru &
Liu (2025) TV-HTE estimator with (2026) feedback extension.

The two packages should produce matching estimates up to MC noise on
the same DGP. Mirror changes both ways.

## What's where

- `src/design.jl` — `unit_design()` and `unit_moments()`. Mirrors R's
  `R/design.R`. Note: `LinearAlgebra.I` is a `UniformScaling`, so the
  PD ridge has to be added element-wise via a `for` loop, not
  `1e-10 .* I` (that's a bug I hit; see commit history).
- `src/simulate.jl` — `simulate_tvhte()`, including the feedback DGP.
- `src/tvhte.jl` — main two-step estimator. Optim.jl + BFGS for QMLE,
  Distributions.jl `MvNormal` for the marginal density. Designs are
  cached per unique cohort under staggered timing.
- `src/feedback.jl` — `fit_feedback()` (pooled OLS via `\`) and
  `simulate_counterfactual()`.
- `test/runtests.jl` — smoke tests at N=300–400.

## Docs

Documenter.jl-based. `docs/make.jl` builds, `.github/workflows/docs.yml`
deploys to `gh-pages` on push. Live at
<https://xiangao.github.io/TVHTE.jl/dev/>.

**Reference page pattern**: every function in `docs/src/reference.md`
gets a `@docs` block followed by an `@example` block (the docstring
shows the API; the example block runs and prints real output). Plain
`# Examples` inside docstrings render as text only — they don't
execute under Documenter. Don't rely on them alone.

## Pages config

GH Pages source = `gh-pages` branch root (one-time `gh api` call;
already done).

## Mirroring R changes

When `tvhte` (R) changes:
1. Make the equivalent change in the matching `src/*.jl` file.
2. Run `julia --project=. test/runtests.jl`.
3. Keep the NamedTuple return shape stable so existing tests and
   vignette examples don't break.
