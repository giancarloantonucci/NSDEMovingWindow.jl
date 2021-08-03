using MovingWindow
using Documenter

DocMeta.setdocmeta!(MovingWindow, :DocTestSetup, :(using MovingWindow); recursive=true)

makedocs(;
    modules=[MovingWindow],
    authors="Giancarlo A. Antonucci",
    repo="https://github.com/antonuccig/MovingWindow.jl/blob/{commit}{path}#{line}",
    sitename="MovingWindow.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://antonuccig.github.io/MovingWindow.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/antonuccig/MovingWindow.jl",
)
