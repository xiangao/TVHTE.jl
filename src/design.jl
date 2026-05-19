# Per-unit linear-Gaussian system implied by Botosaru-Liu (2025), as
# in the R package's R/design.R. Solving the AR(1) in calendar time:
#
#   Y_i = rho_Y (A r) Y_{i0}
#       + (A 1_T)   * alpha_i
#       + (A M c)   * delta_{i0}
#       + (A M W)   * eps_i
#       + A         * U_i
#
# A: (I - rho_Y L)^{-1}, lower triangular with A[t, s] = rho_Y^(t-s).
# M: T x (J+1) indicator placing delta_{i,j} at calendar time t0+j.
# c: (J+1) vector (1, rho_delta, ..., rho_delta^J).
# W: (J+1) x J lower triangular for eps innovations.

struct UnitDesign
    A::Matrix{Float64}
    M::Matrix{Float64}
    c::Vector{Float64}
    W::Matrix{Float64}
    r::Vector{Float64}
    AMc::Vector{Float64}
    A1::Vector{Float64}
    AMW::Matrix{Float64}
end

function unit_design(T::Int, t0, J::Int, rho_Y::Float64, rho_delta::Float64)
    A = zeros(T, T)
    for t in 1:T, s in 1:t
        A[t, s] = rho_Y ^ (t - s)
    end

    M = zeros(T, J + 1)
    for j in 0:J
        t_idx = t0 + j
        if 1 <= t_idx <= T
            M[Int(t_idx), j + 1] = 1.0
        end
    end

    c = [rho_delta ^ j for j in 0:J]

    W = zeros(J + 1, J)
    for j in 1:J, k in 1:j
        W[j + 1, k] = rho_delta ^ (j - k)
    end

    r = zeros(T); r[1] = 1.0

    AMc = vec(A * M * c)
    A1  = vec(A * ones(T))
    AMW = A * M * W

    UnitDesign(A, M, c, W, r, AMc, A1, AMW)
end

# Per-unit marginal mean and covariance of Y_i given Y_{i0}, with optional
# pre-computed X_i β contribution to the mean.
function unit_moments(ds::UnitDesign, Y0::Float64,
                      rho_Y::Float64,
                      sigma_U2::Float64, sigma_eps2::Float64,
                      mu_alpha::Float64, mu_delta0::Float64,
                      sigma_alpha2::Float64, sigma_delta0_2::Float64,
                      cov_alpha_delta::Float64;
                      Xi_beta::Union{Nothing,Vector{Float64}} = nothing)
    mu_Y = rho_Y .* (ds.A * ds.r) .* Y0 .+
           ds.A1  .* mu_alpha .+
           ds.AMc .* mu_delta0
    if Xi_beta !== nothing
        mu_Y .+= Xi_beta
    end

    Sigma_Y = ds.A1 * ds.A1' .* sigma_alpha2 .+
              ds.AMc * ds.AMc' .* sigma_delta0_2 .+
              (ds.A1 * ds.AMc' .+ ds.AMc * ds.A1') .* cov_alpha_delta .+
              ds.AMW * ds.AMW' .* sigma_eps2 .+
              ds.A * ds.A' .* sigma_U2
    Sigma_Y .= (Sigma_Y .+ Sigma_Y') ./ 2
    for i in 1:size(Sigma_Y, 1); Sigma_Y[i, i] += 1e-10; end

    (mean = vec(mu_Y), cov = Symmetric(Sigma_Y))
end
