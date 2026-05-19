"""
    tvhte(Y, Y0; t0, J, X=nothing, compute_se=true)

Two-step QMLE + Tweedie/Gaussian empirical-Bayes estimator for the
Botosaru-Liu (2025) TV-HTE model. See the R companion for details.

# Returns
NamedTuple with `theta`, `prior`, `beta`, `se`, `loglik`,
`lambda_hat`, `delta_path`, `t0`, `J`, `N`, `T`, `K`, `converged`.
"""
function tvhte(Y::Matrix{Float64}, Y0::AbstractVector;
               t0, J::Int, X::Union{Nothing,Array{Float64,3}} = nothing,
               compute_se::Bool = true, maxiter::Int = 500)
    N, T = size(Y)
    length(Y0) == N || throw(ArgumentError("Y0 length mismatch"))
    J >= 0 || throw(ArgumentError("J must be >= 0"))

    t0_vec = isa(t0, AbstractVector) ? Float64.(t0) : fill(Float64(t0), N)
    length(t0_vec) == N || throw(ArgumentError("t0 must be scalar or length N"))

    K = X === nothing ? 0 : size(X, 3)
    if X !== nothing
        size(X, 1) == N && size(X, 2) == T || throw(ArgumentError("X must be N x T x K"))
    end

    # Parameter packing (real-valued; tanh/exp on natural-scale params)
    # Layout: [atanh(rho_Y), atanh(rho_delta), log(sigma_U2), log(sigma_eps2),
    #          mu_alpha, mu_delta0, log(sigma_alpha2), log(sigma_delta0_2),
    #          atanh(cor), beta1, ..., beta_K]
    var_Y = var(vec(Y))
    p0 = [atanh(0.3); atanh(0.3); log(var_Y); log(0.25 * var_Y);
          mean(Y); 0.0; log(0.25 * var_Y); log(0.25 * var_Y);
          0.0; zeros(K)]

    unpack(p) = (rho_Y = tanh(p[1]),
                 rho_delta = tanh(p[2]),
                 sigma_U2 = exp(p[3]),
                 sigma_eps2 = exp(p[4]),
                 mu_alpha = p[5], mu_delta0 = p[6],
                 sigma_alpha2 = exp(p[7]),
                 sigma_delta0_2 = exp(p[8]),
                 cov_alpha_delta = tanh(p[9]) * sqrt(exp(p[7]) * exp(p[8])),
                 beta = K == 0 ? Float64[] : p[10:(9 + K)])

    unique_t0 = unique(t0_vec)
    build_ds(rho_Y, rho_delta) = Dict(t => unit_design(T, t, J, rho_Y, rho_delta)
                                       for t in unique_t0)

    function unit_Xbeta(ds_dict, beta_vec)
        K == 0 && return nothing
        A = first(values(ds_dict)).A
        [vec(A * (reshape(X[i, :, :], T, K) * beta_vec)) for i in 1:N]
    end

    function nll(p)
        par = unpack(p)
        ds_dict = build_ds(par.rho_Y, par.rho_delta)
        XB = unit_Xbeta(ds_dict, par.beta)
        ll = 0.0
        for i in 1:N
            ds_i = ds_dict[t0_vec[i]]
            Xi_beta = XB === nothing ? nothing : XB[i]
            mom = unit_moments(ds_i, Float64(Y0[i]),
                               par.rho_Y, par.sigma_U2, par.sigma_eps2,
                               par.mu_alpha, par.mu_delta0,
                               par.sigma_alpha2, par.sigma_delta0_2,
                               par.cov_alpha_delta;
                               Xi_beta = Xi_beta)
            d = try
                MvNormal(mom.mean, mom.cov)
            catch
                return prevfloat(Inf)
            end
            ll_i = logpdf(d, vec(Y[i, :]))
            if !isfinite(ll_i)
                return prevfloat(Inf)
            end
            ll += ll_i
        end
        -ll
    end

    opt_res = optimize(nll, p0, BFGS(), Options(iterations = maxiter))
    phat = minimizer(opt_res)
    par = unpack(phat)
    loglik = -opt_res.minimum

    # Step 2: posterior mean of lambda_i (Gaussian conjugacy)
    ds_dict = build_ds(par.rho_Y, par.rho_delta)
    XB = unit_Xbeta(ds_dict, par.beta)

    mu_lambda = [par.mu_alpha, par.mu_delta0]
    lambda_hat = zeros(N, 2)
    for i in 1:N
        ds_i = ds_dict[t0_vec[i]]
        Xi_beta = XB === nothing ? nothing : XB[i]
        mom = unit_moments(ds_i, Float64(Y0[i]),
                           par.rho_Y, par.sigma_U2, par.sigma_eps2,
                           par.mu_alpha, par.mu_delta0,
                           par.sigma_alpha2, par.sigma_delta0_2,
                           par.cov_alpha_delta; Xi_beta = Xi_beta)
        Sigma_lY = vcat(
            (ds_i.A1 .* par.sigma_alpha2 .+ ds_i.AMc .* par.cov_alpha_delta)',
            (ds_i.A1 .* par.cov_alpha_delta .+ ds_i.AMc .* par.sigma_delta0_2)'
        )
        inv_S_Y = inv(Matrix(mom.cov))
        lambda_hat[i, :] = mu_lambda .+ Sigma_lY * inv_S_Y * (vec(Y[i, :]) .- mom.mean)
    end

    # Posterior event-time trajectory delta_{i, j}
    delta_path = zeros(N, J + 1)
    delta_path[:, 1] .= lambda_hat[:, 2]
    if J >= 1
        for i in 1:N
            ds_i = ds_dict[t0_vec[i]]
            Xi_beta = XB === nothing ? nothing : XB[i]
            mom = unit_moments(ds_i, Float64(Y0[i]),
                               par.rho_Y, par.sigma_U2, par.sigma_eps2,
                               par.mu_alpha, par.mu_delta0,
                               par.sigma_alpha2, par.sigma_delta0_2,
                               par.cov_alpha_delta; Xi_beta = Xi_beta)
            Sigma_epsY = transpose(ds_i.AMW) .* par.sigma_eps2
            inv_S_Y = inv(Matrix(mom.cov))
            eps_hat = Sigma_epsY * inv_S_Y * (vec(Y[i, :]) .- mom.mean)
            for j in 1:J
                delta_path[i, j + 1] = par.rho_delta ^ j * lambda_hat[i, 2] +
                                       sum(par.rho_delta .^ ((j - 1):-1:0) .* eps_hat[1:j])
            end
        end
    end

    (theta = (rho_Y = par.rho_Y, rho_delta = par.rho_delta,
              sigma_U2 = par.sigma_U2, sigma_eps2 = par.sigma_eps2),
     prior = (mu_alpha = par.mu_alpha, mu_delta0 = par.mu_delta0,
              sigma_alpha2 = par.sigma_alpha2,
              sigma_delta0_2 = par.sigma_delta0_2,
              cov_alpha_delta = par.cov_alpha_delta),
     beta = par.beta,
     loglik = loglik,
     converged = converged(opt_res),
     lambda_hat = lambda_hat, delta_path = delta_path,
     t0 = t0_vec, J = J, N = N, T = T, K = K)
end
