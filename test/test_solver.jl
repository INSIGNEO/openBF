using Base.Test
using openBF

@testset "solver.jl" begin

    project_name = "test"
    inputs = openBF.loadSimulationFiles(project_name)
    vessels, edges = openBF.buildArterialNetwork(inputs[2], inputs[3][1], inputs[4])

    Ccfl = inputs[1]["Ccfl"]
    dt = openBF.calculateDeltaT(vessels, Ccfl)
    @test isapprox(dt, 1e-4, atol=1e-4); println("calculateDeltaT - [OK]")

    @test_nowarn openBF.solveModel(vessels, inputs[3], edges, inputs[4], dt, 0.0); println("solveModel - [OK]")

    v = vessels[1]
    @test_nowarn openBF.computeLimiter(v, v.vA, v.invDx, v.dU, v.slopesA); println("computeLimiter - [OK]")


end
