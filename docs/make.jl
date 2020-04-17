using Documenter, GridapODEs

makedocs(;
    modules=[GridapODEs],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/santiagobadia/GridapODEs.jl/blob/{commit}{path}#L{line}",
    sitename="GridapODEs.jl",
    authors="Santiago Badia",
    assets=String[],
)

deploydocs(;
    repo="github.com/santiagobadia/GridapODEs.jl",
)
