using Base.Test
using openBF

@testset "boundary_conditions.jl" begin

    project_name = "test"
    inputs = openBF.loadSimulationFiles(project_name)
    heart = inputs[3][1]
    vessels, edge_list = openBF.buildArterialNetwork(inputs[2], heart, inputs[4])

    Pref = vessels[1].P[1]
    Qref = vessels[1].Q[1]
    Ccfl = inputs[1]["Ccfl"]
    dt = openBF.calculateDeltaT(vessels, Ccfl)

    @test heart.inlet_type == "Q"
    @inferred openBF.setInletBC(0.0, dt, vessels[1], heart)
    @test vessels[1].Q[1] != Qref
    @test vessels[1].P[1] != Pref

    @inferred openBF.inletCompatibility(dt, vessels[1], heart)
    @test isapprox(vessels[1].A[1], 1.8e-4, atol=1e-5)
    @test isapprox(vessels[1].P[1], -1.6e-9, atol=1e-9)

    heart.inlet_type = "P"
    openBF.setInletBC(0.0, dt, vessels[1], heart)
    @test vessels[1].Q[1] != Qref
    @test vessels[1].P[1] != Pref

    @inferred openBF.inletCompatibility(dt, vessels[1], heart)
    @test isapprox(vessels[1].A[1], 1.8e-4, atol=1e-5)
    @test isapprox(vessels[1].Q[1], -5.2e-7, atol=1e-7)

    @test isapprox(openBF.inputFromData(0.0, heart), -5.24e-07, atol=1e-7)

    @test_nowarn openBF.setOutletBC(dt, vessels[3])
    @test isapprox(vessels[3].A[end], 9.5e-5, atol=1e-5)
    @test vessels[3].u[end] == 0.0

    M = vessels[1].M
    @inferred openBF.updateGhostCells(vessels[1])
    @test vessels[1].U00A == vessels[1].A[1]
    @test vessels[1].U00Q == vessels[1].Q[1]
    @test vessels[1].U01A == vessels[1].A[2]
    @test vessels[1].U01Q == vessels[1].Q[2]
    @test vessels[1].UM1A == vessels[1].A[M]
    @test vessels[1].UM1Q == vessels[1].Q[M]
    @test vessels[1].UM2A == vessels[1].A[M-1]
    @test vessels[1].UM2Q == vessels[1].Q[M-1]

end
