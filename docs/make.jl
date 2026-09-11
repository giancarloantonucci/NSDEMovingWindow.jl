using Documenter
using NSDEMovingWindow

PAGES = [
    "Home" => "index.md",
    "Strategies" => "strategies.md",
    "API" => "api.md"
]

makedocs(;
    sitename = "NSDEMovingWindow.jl",
    format = Documenter.HTML(),
    modules = [NSDEMovingWindow],
    pages = PAGES,
    checkdocs = :exports, # every export must carry a docstring, or the build fails
    authors = "Giancarlo A. Antonucci <giancarlo.antonucci@icloud.com>"
)

deploydocs(;
    repo = "github.com/giancarloantonucci/NSDEMovingWindow.jl.git"
)
