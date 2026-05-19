using Pkg
Pkg.develop(PackageSpec(path = joinpath(@__DIR__, "..")))
Pkg.instantiate()

using Documenter
using TVHTE

makedocs(
    sitename = "TVHTE.jl",
    modules  = [TVHTE],
    format   = Documenter.HTML(
        edit_link = nothing,
        repolink  = "https://github.com/xiangao/TVHTE.jl",
    ),
    pages = [
        "Home"      => "index.md",
        "Vignettes" => ["Getting Started" => "vignettes/01_getting_started.md"],
        "Reference" => "reference.md",
    ],
    warnonly  = true,
    checkdocs = :none,
    remotes   = nothing,
)

deploydocs(
    repo         = "github.com/xiangao/TVHTE.jl.git",
    devbranch    = "master",
    push_preview = false,
)
