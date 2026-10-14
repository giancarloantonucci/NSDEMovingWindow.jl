# Run collectively with at least two ranks and the NSDE packages available.
using NSDEMovingWindow, NSDERungeKutta, Test
using NSDEBase: NSDEBase
const MPI = NSDETimeParallel.MPI
MPI.Init()
const rank = MPI.Comm_rank(MPI.COMM_WORLD)
const N = MPI.Comm_size(MPI.COMM_WORLD)
@testset "MPI windows inherit the seeded initial state" begin
    problem = Logistic(0.3, (0, 3))
    parareal = Parareal(RK4(h=0.01), Euler(h=0.1); N, criterion=NoStop())
    mowi = MoWi(parareal; τ=1.0, Δτ=0.5)
    for backend in (MPIBackend(), PipelinedMPIBackend())
        cache = NSDEMovingWindow.MoWiCache(problem, mowi)
        result = MoWiSolution(problem, mowi)
        mowi(cache, result, problem; backend)
        serial = solve(problem, mowi)
        @test result.status == :converged
        @test length(result) == length(serial)
        @test [boundaryvalues(w.lastiterate) for w in result.windows] == [boundaryvalues(w.lastiterate) for w in serial.windows]
        @test result[2][1].u[1] != problem.u0
        @test result[2][1].u[1] == result[1](0.5)
        q = remake(problem; u0=[0.7])
        fill!(result.restarts, 99)
        mowi(cache, result, q; backend)
        fresh = solve(q, mowi)
        @test result.restarts == fresh.restarts
        @test [boundaryvalues(w.lastiterate) for w in result.windows] == [boundaryvalues(w.lastiterate) for w in fresh.windows]
    end
end
@testset "MPI nonfinite window stops collectively" begin
    p = Dahlquist([1.0], (0.0, 32.0); λ=-1000.0)
    parareal = Parareal(RK4(h=0.001), Euler(h=0.02); N)
    mowi = MoWi(parareal; τ=16.0, Δτ=8.0)
    result = @test_logs (:warn, r"stopped at a nonfinite") solve(p, mowi; backend=MPIBackend())
    @test result.status == :nonfinite
    @test length(result) == 1
    @test isempty(result[1].errors)
    @test_throws ArgumentError ensemblemean(result)
end

@testset "Decimal endpoint policy agrees on every rank" begin
    for backend in (MPIBackend(), PipelinedMPIBackend())
        p=Logistic(0.3, (0.0, 1.0))
        par=Parareal(RK4(h=0.005), Euler(h=0.01); N=4, criterion=NoStop())
        sol=solve(p, MoWi(par; τ=0.1, Δτ=0.1); backend)
        @test sol.status == :converged
        @test length(sol) == 10
        @test NSDEMovingWindow.windowspan(sol[end])[2] == 1.0
        @test !any(sol.clipped)
    end
end
MPI.Finalize()
