using Test
using TVHTE
using LinearAlgebra
using Statistics
using Random

@testset "simulate + tvhte basic" begin
    sim = simulate_tvhte(N = 400, T = 6, t0 = 3, J = 3,
                         rho_Y = 0.5, rho_delta = 0.7,
                         sigma_alpha = 0.5, sigma_delta0 = 0.4,
                         seed = 1)
    fit = tvhte(sim.Y, sim.Y0; t0 = sim.t0, J = sim.J,
                compute_se = false)
    @test fit.converged
    @test abs(fit.theta.rho_Y     - 0.5) < 0.15
    @test abs(fit.theta.rho_delta - 0.7) < 0.20
end

@testset "covariates" begin
    sim = simulate_tvhte(N = 400, T = 6, t0 = 3, J = 3,
                         beta = [0.7, -0.4], seed = 3)
    fit = tvhte(sim.Y, sim.Y0; t0 = sim.t0, J = sim.J, X = sim.X,
                compute_se = false)
    @test length(fit.beta) == 2
    @test maximum(abs.(fit.beta .- [0.7, -0.4])) < 0.20
end

@testset "staggered cohorts incl Inf" begin
    Random.seed!(20)
    N = 400
    cohorts = rand([3.0, 4.0, 5.0, Inf], N)
    sim = simulate_tvhte(N = N, T = 7, t0 = cohorts, J = 2, seed = 21)
    fit = tvhte(sim.Y, sim.Y0; t0 = sim.t0, J = sim.J,
                compute_se = false)
    @test fit.converged
    @test abs(fit.theta.rho_Y - 0.5) < 0.25
end

@testset "feedback + counterfactual" begin
    sim = simulate_tvhte(N = 400, T = 6, t0 = 3, J = 2,
                         beta = 0.4,
                         feedback_gamma = [0.2, 0.3, 0.5, 0.4],
                         seed = 51)
    fb = fit_feedback(sim.Y, sim.Y0, sim.X, sim.X0)
    @test abs(fb.coef.intercept - 0.2) < 0.15
    @test abs(fb.coef.gamma_Y   - 0.3) < 0.15
    @test abs(fb.coef.gamma_X   - 0.5) < 0.15

    fit = tvhte(sim.Y, sim.Y0; t0 = sim.t0, J = sim.J, X = sim.X,
                compute_se = false)
    cf = simulate_counterfactual(fit, fb, 4; N_star = 200, seed = 99)
    @test size(cf.Y) == (200, 6)
    @test size(cf.X) == (200, 6, 1)
    @test all(isfinite, cf.Y)
end
