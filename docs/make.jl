using Documenter
using MovingWindow

PAGES = ["Home" => "index.md"]

makedocs(;
    sitename = "MovingWindow.jl",
    format = Documenter.HTML(),
    modules = [MovingWindow],
    pages = PAGES,
    authors = "Giancarlo A. Antonucci <giancarlo.antonucci@icloud.com>"
)

deploydocs(;
    repo = "https://github.com/antonuccig/MovingWindow.jl"
)
