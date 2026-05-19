"""
    fit_feedback(Y, Y0, X, X0)

Pooled OLS for the homogeneous covariate-feedback process
(Botosaru-Liu 2026). Mirrors `fit_feedback()` in the R package.
"""
function fit_feedback(Y::Matrix{Float64}, Y0::AbstractVector,
                      X::Array{Float64,3}, X0::AbstractVector)
    size(X, 3) == 1 || throw(ArgumentError("Phase 5 scope is K = 1 covariate"))
    N, T = size(Y)
    length(Y0) == N || throw(ArgumentError("Y0 length mismatch"))
    length(X0) == N || throw(ArgumentError("X0 length mismatch"))

    # Stack X_it, Y_{i,t-1}, X_{i,t-1} across all (i, t)
    X_flat = vec(X[:, :, 1])
    Y_lag  = vec(hcat(Y0, Y[:, 1:T-1]))
    X_lag  = vec(hcat(X0, X[:, 1:T-1, 1]))

    XX = hcat(ones(length(X_flat)), Y_lag, X_lag)
    coef = (XX' * XX) \ (XX' * X_flat)
    resid = X_flat .- XX * coef
    sigma_eta = sqrt(mean(resid .^ 2))

    (coef = (intercept = coef[1], gamma_Y = coef[2], gamma_X = coef[3]),
     sigma_eta = sigma_eta)
end

"""
    simulate_counterfactual(fit, fit_fb, t0_star; N_star=500, Y0_star=nothing, X0_star=nothing, seed=nothing)

Algorithm 1 of Botosaru-Liu (2026): joint counterfactual `(Y*, X*)`
under alternative treatment timing, using the estimated structural and
feedback models.
"""
function simulate_counterfactual(fit, fit_fb, t0_star;
                                 N_star::Int = 500,
                                 Y0_star::Union{Nothing,AbstractVector} = nothing,
                                 X0_star::Union{Nothing,AbstractVector} = nothing,
                                 seed::Union{Nothing,Int} = nothing)
    seed === nothing || Random.seed!(seed)
    length(fit.beta) == 1 || throw(ArgumentError("Phase 5 scope is K = 1 covariate"))

    T = fit.T; J = fit.J
    t0_vec = isa(t0_star, AbstractVector) ? Float64.(t0_star) :
             fill(Float64(t0_star), N_star)
    length(t0_vec) == N_star || throw(ArgumentError("t0_star length mismatch"))

    Y0v = Y0_star === nothing ? randn(N_star) : Float64.(Y0_star)
    X0v = X0_star === nothing ? randn(N_star) : Float64.(X0_star)

    Sigma = [fit.prior.sigma_alpha2     fit.prior.cov_alpha_delta;
             fit.prior.cov_alpha_delta  fit.prior.sigma_delta0_2]
    Sigma_pd = copy(Sigma)
    for i in 1:size(Sigma_pd, 1); Sigma_pd[i, i] += 1e-10; end
    L = cholesky(Symmetric(Sigma_pd)).L
    Z = randn(N_star, 2)
    lambda = Z * L' .+ [fit.prior.mu_alpha fit.prior.mu_delta0]
    alpha  = lambda[:, 1]
    delta0 = lambda[:, 2]

    sigma_eps = sqrt(fit.theta.sigma_eps2)
    delta = zeros(N_star, J + 1)
    delta[:, 1] .= delta0
    if J >= 1
        eps = randn(N_star, J) .* sigma_eps
        for j in 1:J
            delta[:, j + 1] .= fit.theta.rho_delta .* delta[:, j] .+ eps[:, j]
        end
    end

    beta = fit.beta[1]
    sigma_U = sqrt(fit.theta.sigma_U2)
    rho_Y = fit.theta.rho_Y
    g0 = fit_fb.coef.intercept
    gY = fit_fb.coef.gamma_Y
    gX = fit_fb.coef.gamma_X
    sigma_eta = fit_fb.sigma_eta

    Y = zeros(N_star, T)
    X = zeros(N_star, T, 1)
    for t in 1:T
        Y_lag = t == 1 ? Y0v : Y[:, t - 1]
        X_lag = t == 1 ? X0v : X[:, t - 1, 1]
        X[:, t, 1] .= g0 .+ gY .* Y_lag .+ gX .* X_lag .+ randn(N_star) .* sigma_eta

        j_i = t .- t0_vec
        in_window = isfinite.(j_i) .& (j_i .>= 0) .& (j_i .<= J)
        trt = zeros(N_star)
        for i in 1:N_star
            if in_window[i]
                trt[i] = delta[i, Int(j_i[i]) + 1]
            end
        end

        Y[:, t] .= rho_Y .* Y_lag .+ alpha .+ trt .+
                   X[:, t, 1] .* beta .+ randn(N_star) .* sigma_U
    end

    (Y = Y, X = X, lambda = lambda, t0 = t0_vec)
end
