module TVHTE

using LinearAlgebra
using Statistics
using Random
using Distributions: MvNormal, Normal, logpdf
using Optim: optimize, BFGS, minimizer, converged, Options

export simulate_tvhte, tvhte, fit_feedback, simulate_counterfactual

include("design.jl")
include("simulate.jl")
include("tvhte.jl")
include("feedback.jl")

end # module
