"""
    simulate_tvhte(; N, T, t0, J, ...)

Simulate from the Botosaru-Liu (2025) TV-HTE model with optional
covariate feedback (Botosaru-Liu 2026). Mirrors the R package's
`simulate_tvhte()`.

# Keywords
- `t0`: scalar (common timing) or length-`N` vector. `Inf` for never-treated.
- `J`: max event time observed in-window.
- `rho_Y`, `rho_delta`: outcome and event-time AR(1) coefficients.
- `sigma_U`, `sigma_eps`: SDs of `U_it` and `eps_ij`.
- `mu_alpha`, `mu_delta0`, `sigma_alpha`, `sigma_delta0`,
  `cor_alpha_delta`: Gaussian prior on `lambda_i`.
- `Y0_mean`, `Y0_sd`: baseline outcome.
- `beta`: optional scalar (or length-K) coefficient on covariate(s).
- `feedback_gamma`: optional length-4 vector `(γ0, γY, γX, σ_η)`
  enabling endogenous AR(1) covariate dynamics (K = 1 only).
- `seed`: RNG seed.
"""
function simulate_tvhte(; N::Int = 500, T::Int = 6,
                          t0 = 3, J::Int = 3,
                          rho_Y::Float64 = 0.5, rho_delta::Float64 = 0.7,
                          sigma_U::Float64 = 1.0, sigma_eps::Float64 = 0.3,
                          mu_alpha::Float64 = 0.0, mu_delta0::Float64 = 1.0,
                          sigma_alpha::Float64 = 0.5, sigma_delta0::Float64 = 0.5,
                          cor_alpha_delta::Float64 = 0.0,
                          Y0_mean::Float64 = 0.0, Y0_sd::Float64 = 1.0,
                          beta::Union{Nothing,Real,AbstractVector} = nothing,
                          feedback_gamma::Union{Nothing,AbstractVector} = nothing,
                          seed::Union{Nothing,Int} = nothing)
    seed === nothing || Random.seed!(seed)
    J >= 0 || throw(ArgumentError("J must be >= 0"))

    t0_vec = isa(t0, AbstractVector) ? Float64.(t0) : fill(Float64(t0), N)
    length(t0_vec) == N || throw(ArgumentError("t0 must be scalar or length N"))
    for v in t0_vec
        if isfinite(v) && (v < 1 || v > T)
            throw(ArgumentError("finite t0 entries must be in 1:T"))
        end
    end

    # lambda ~ N(mu, Sigma_lambda)
    Sigma_lambda = [sigma_alpha^2                            cor_alpha_delta*sigma_alpha*sigma_delta0;
                    cor_alpha_delta*sigma_alpha*sigma_delta0 sigma_delta0^2]
    L = cholesky(Sigma_lambda).L
    Z = randn(N, 2)
    lambda = Z * L' .+ [mu_alpha mu_delta0]
    alpha  = lambda[:, 1]
    delta0 = lambda[:, 2]

    # delta_{i,j} via AR(1) starting from delta_{i,0}
    eps = randn(N, J) .* sigma_eps
    delta = zeros(N, J + 1)
    delta[:, 1] .= delta0
    for j in 1:J
        delta[:, j + 1] .= rho_delta .* delta[:, j] .+ eps[:, j]
    end

    Y0 = Y0_mean .+ Y0_sd .* randn(N)

    # Covariate setup
    beta_vec = beta === nothing ? Float64[] :
               (isa(beta, AbstractVector) ? Float64.(beta) : [Float64(beta)])
    K = length(beta_vec)
    X  = K == 0 ? nothing : zeros(N, T, K)
    X0 = feedback_gamma === nothing ? nothing : randn(N)
    if K > 0 && feedback_gamma === nothing
        X .= randn(N, T, K)
    end
    if feedback_gamma !== nothing && K != 1
        throw(ArgumentError("feedback_gamma requires K == 1"))
    end

    Y = zeros(N, T)
    for t in 1:T
        Y_lag = t == 1 ? Y0 : Y[:, t - 1]

        if feedback_gamma !== nothing
            X_lag = t == 1 ? X0 : X[:, t - 1, 1]
            X[:, t, 1] .= feedback_gamma[1] .+
                          feedback_gamma[2] .* Y_lag .+
                          feedback_gamma[3] .* X_lag .+
                          randn(N) .* feedback_gamma[4]
        end

        j_i = t .- t0_vec
        in_window = isfinite.(j_i) .& (j_i .>= 0) .& (j_i .<= J)
        trt = zeros(N)
        for i in 1:N
            if in_window[i]
                trt[i] = delta[i, Int(j_i[i]) + 1]
            end
        end

        x_eff = zeros(N)
        if K > 0
            for i in 1:N
                x_eff[i] = dot(view(X, i, t, :), beta_vec)
            end
        end

        Y[:, t] .= rho_Y .* Y_lag .+ alpha .+ trt .+ x_eff .+ randn(N) .* sigma_U
    end

    (Y = Y, Y0 = Y0, X = X, X0 = X0, t0 = t0_vec, J = J,
     lambda = lambda, delta = delta,
     params = (rho_Y = rho_Y, rho_delta = rho_delta,
               sigma_U = sigma_U, sigma_eps = sigma_eps,
               mu_alpha = mu_alpha, mu_delta0 = mu_delta0,
               sigma_alpha = sigma_alpha, sigma_delta0 = sigma_delta0,
               cor_alpha_delta = cor_alpha_delta,
               beta = beta_vec,
               feedback_gamma = feedback_gamma))
end
