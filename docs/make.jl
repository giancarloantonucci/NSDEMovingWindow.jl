using Documenter
using NSDEMovingWindow

PAGES = ["Home" => "index.md"]

makedocs(;
    sitename = "NSDEMovingWindow.jl",
    format = Documenter.HTML(),
    modules = [NSDEMovingWindow],
    pages = PAGES,
    authors = "Giancarlo A. Antonucci <giancarlo.antonucci@icloud.com>"
)

deploydocs(;
    repo = "https://github.com/giancarloantonucci/NSDEMovingWindow.jl"
)
