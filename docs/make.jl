using Documenter, GridapTimeStepper

makedocs(;
    modules=[GridapTimeStepper],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/santiagobadia/GridapTimeStepper.jl/blob/{commit}{path}#L{line}",
    sitename="GridapTimeStepper.jl",
    authors="Santiago Badia",
    assets=String[],
)

deploydocs(;
    repo="github.com/santiagobadia/GridapTimeStepper.jl",
)
