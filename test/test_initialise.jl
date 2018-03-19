using Base.Test
using openBF

@testset "initialise.jl" begin

    project_name = "test"
    first_inlet = "$project_name\_inlet.dat"
    inlets = [first_inlet]
    number_of_inlets = 2

    @test_nowarn openBF.checkInputFiles(project_name)

    inlets = openBF.checkInletFiles(project_name, number_of_inlets, inlets)
    @test length(inlets) == number_of_inlets

    @test_nowarn openBF.copyInputFilesToResultsFolder(project_name, inlets)
    @test isfile("$project_name.csv")
    @test isfile("$project_name\_constants.yml")
    @test isfile("$project_name\_inlet.dat")
    @test isfile("$project_name\_2_inlet.dat")

    model, model_header = openBF.readModelData("$project_name.csv")
    @test typeof(model) == Array{Any,2}
    @test model_header == ["name" "sn" "tn" "in" "L" "M" "Rp" "Rd" "E" "Pext" "R1" "R2" "C"]

    constants = openBF.loadConstants("$project_name\_constants.yml")
    @test typeof(constants) == Dict{Any,Any} #constants

    @test openBF.checkConstants(constants) == constants

    blood = openBF.buildBlood(constants)
    @test typeof(blood) == Blood

    inlet_data = openBF.loadInletData("$project_name\_inlet.dat")
    @test typeof(inlet_data) == Array{Float64,2}

    heart = openBF.buildHeart(constants, inlet_data, 1)
    @test typeof(heart) == Heart

    @inferred openBF.buildHearts(constants, inlets)

    dx, invDx, halfDx = openBF.meshVessel(1.0, 10)
    @test dx == 0.001
    @test invDx == 1000.0
    @test halfDx == 0.0005

    A0 = ones(Float64, 2)
    gamma = copy(A0)
    BCout = ["", "", "", ""]
    BCout, outlet = openBF.detectCapillaries(BCout, blood, A0, gamma)
    @test outlet == "none"
    @test BCout[1] == 0.0
    @test BCout[2] == 0.0
    @test BCout[3] == 0.0
    @test BCout[4] == 0.0

    BCout = [1.0, "", "", ""]
    BCout, outlet = openBF.detectCapillaries(BCout, blood, A0, gamma)
    @test outlet == "reflection"
    @test BCout[1] == 1.0
    @test BCout[2] == 0.0
    @test BCout[3] == 0.0
    @test BCout[4] == 0.0

    BCout = ["", "", 1.0, 1.0]
    BCout, outlet = openBF.detectCapillaries(BCout, blood, A0, gamma)
    @test outlet == "wk3"
    @test BCout[1] == 0.0
    @test isapprox(BCout[2],  1298, atol=1.0)
    @test isapprox(BCout[3], -1297, atol=1.0)
    @test BCout[4] == 1.0

    BCout = ["", 1.0, 1.0, 1.0]
    BCout, outlet = openBF.detectCapillaries(BCout, blood, A0, gamma)
    @test outlet == "wk3"
    @test BCout[1] == 0.0
    @test BCout[2] == 1.0
    @test BCout[3] == 1.0
    @test BCout[3] == 1.0

    m_row = model[1,:]
    @test_nowarn openBF.parseModelRow(m_row)
    vessel_name, sn, tn, rn, L, M, Rp, Rd, E, Pext, BCout =  openBF.parseModelRow(m_row)

    vessel = openBF.buildVessel(1, m_row, heart, blood)
    @test typeof(vessel) == Vessel

    vessels, edge_list = openBF.buildArterialNetwork(model, heart, blood)
    @test typeof(vessels) == Array{Vessel,1}
    @test typeof(edge_list) == Array{Int64,2}

    cd("..")
    rm("$project_name\_results", recursive=true)

    inputs = openBF.loadSimulationFiles(project_name)
    @test typeof(inputs[1]) == Dict{Any,Any} #constants
    @test typeof(inputs[2]) == Array{Any,2} #model
    @test typeof(inputs[3]) == Array{Any,1} #hearts
    @test typeof(inputs[4]) == Blood
    @test typeof(inputs[5]) == Float64 #total_time

    cd("..")
    rm("$project_name\_results", recursive=true)
end
