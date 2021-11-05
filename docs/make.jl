using FEMSolids
using Documenter

DocMeta.setdocmeta!(FEMSolids, :DocTestSetup, :(using FEMSolids); recursive=true)

makedocs(;
    modules=[FEMSolids],
    authors="kimauth <auth@chalmers.se> and contributors",
    repo="https://github.com/kimauth/FEMSolids.jl/blob/{commit}{path}#{line}",
    sitename="FEMSolids.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kimauth.github.io/FEMSolids.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kimauth/FEMSolids.jl.git",
    devbranch="main",
    push_preview=true,
)
