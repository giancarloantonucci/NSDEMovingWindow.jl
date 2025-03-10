using Documenter
using DocumenterInterLinks
using NSDEBase, NSDETimeParallel, NSDEMovingWindow

PAGES = ["Home" => "index.md"]

links = InterLinks(
    "NSDEBase" => (
        "https://giancarloantonucci.github.io/NSDEBase.jl/dev/",
        "https://giancarloantonucci.github.io/NSDEBase.jl/dev/objects.inv"
    ),
    "NSDETimeParallel" => (
        "https://giancarloantonucci.github.io/NSDETimeParallel.jl/dev/",
        "https://giancarloantonucci.github.io/NSDETimeParallel.jl/dev/objects.inv"
    )
)

makedocs(;
    sitename = "NSDEMovingWindow.jl",
    format = Documenter.HTML(),
    modules = [NSDEMovingWindow],
    pages = PAGES,
    authors = "Giancarlo A. Antonucci <giancarlo.antonucci@icloud.com>",
    plugins = [links],
)

deploydocs(;
    repo = "https://github.com/giancarloantonucci/NSDEMovingWindow.jl"
)
