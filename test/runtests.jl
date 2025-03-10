using NSDEMovingWindow
using Test

using Pkg
for (uuid, pkg) in Pkg.dependencies()
    if pkg.name == "NSDEBase"
        @info "Testing NSDERungeKutta with NSDEBase version: $(pkg.version) from $(pkg.source)"
    end
end

@testset "NSDEMovingWindow.jl" begin
    # Write your tests here.
end
