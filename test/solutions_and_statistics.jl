@testset "flatten" begin
    mowi = MoWi(good_parareal(); τ=1.0, Δτ=0.5)
    solution = solve(problem, mowi)
    u, t = flatten(solution)
    @test length(u) == length(t)
    @test length(u) ==
        sum(length(NSDETimeParallel.flatten(solution[m])[1]) for m in 1:length(solution))
    @test !issorted(t) || length(solution) == 1 # window seams rewind time by design
end

@testset "ensemble statistics" begin
    # On the fixed point u ≡ 1 every window is the constant solution: the
    # ensemble mean must be exactly 1, the spread exactly 0.
    fixedpoint = Logistic(1.0, (0.0, 3.0))
    solution = solve(fixedpoint, MoWi(good_parareal(); τ=1.0, Δτ=0.5))
    @test ensemblemean(solution; samples=7) ≈ ones(7) atol = 1e-10
    @test ensemblevariance(solution; samples=7) ≈ zeros(7) atol = 1e-20
    @test ensemblesem(solution; samples=7) ≈ zeros(7) atol = 1e-20
    # General shape and the σ/√M relation on a non-trivial run:
    solution = solve(problem, MoWi(good_parareal(); τ=1.0, Δτ=0.5))
    σ² = ensemblevariance(solution, u -> u[1]; samples=11)
    @test all(≥(0), σ²)
    @test ensemblesem(solution, u -> u[1]; samples=11) ≈ sqrt.(σ² ./ (length(solution) - 1))
end

@testset "evaluator edges (diagnostic use)" begin
    solution = solve(problem, MoWi(good_parareal(); τ=1.0, Δτ=0.5))
    @test solution(-5.0) == solution[1](-5.0)                  # clamped below
    @test solution(99.0) == solution[end](99.0)                # clamped above
    @test solution(0.75) == solution[2](0.75)                  # LATEST covering window wins
    @test norm(solution(1.7) - finetruth(1.7)) < 1e-6
end

@testset "container types" begin
    solution = solve(problem, MoWi(good_parareal(); τ=1.0, Δτ=0.5))
    @test eltype(solution.restarts) == Int
    for m in 1:length(solution)
        @test isconcretetype(eltype(solution[m].lastiterate.chunks)) # the hot seam stays concrete
    end
end

@testset "RecipesBase recipes (headless)" begin
    solution = solve(problem, MoWi(good_parareal(); τ=1.0, Δτ=0.5))
    for obj in (solution, NSDEBase._PhasePlot(solution), NSDEBase._Convergence(solution))
        @test RecipesBase.apply_recipe(Dict{Symbol,Any}(), obj) isa Vector{RecipesBase.RecipeData}
    end
    # The MoWi convergence recipe's 2-element-label path (K_m and R_m
    # series), exercised on purpose:
    @test RecipesBase.apply_recipe(
        Dict{Symbol,Any}(:label => ["K", "R"]), NSDEBase._Convergence(solution)
    ) isa Vector{RecipesBase.RecipeData}
end

@testset "MoWi plot labels" begin
    solution = solve(problem, MoWi(good_parareal(); τ=1.0, Δτ=0.5))
    for (given, expected) in (
        ("MoWi", ["MoWi", ""]),
        ("x", ["x", ""]),
        (:m, ["m", ""]),
        (["a" "b"], ["a", "b"]),
        (["a"], ["a", ""]),
    )
        series = RecipesBase.apply_recipe(
            Dict{Symbol,Any}(:label => given), NSDEBase._Convergence(solution)
        )
        @test length(series) == 2
        @test [s.plotattributes[:label] for s in series] == expected
    end
    phase = RecipesBase.apply_recipe(
        Dict{Symbol,Any}(:label => "MoWi"), NSDEBase._PhasePlot(solution)
    )
    @test phase[1].plotattributes[:label] == "MoWi"
    @test all(s -> s.plotattributes[:label] == "", phase[2:end])
end
