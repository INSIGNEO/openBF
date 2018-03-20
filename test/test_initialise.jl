using Base.Test
using openBF

@testset "initialise.jl" begin

    # project_name = "test"
    # first_inlet = "$project_name\_inlet.dat"
    # inlets = [first_inlet]
    # number_of_inlets = 2
    #
    # @test_nowarn openBF.checkInputFiles(project_name); println("checkInputFiles - [OK]")
    #
    # inlets = openBF.checkInletFiles(project_name, number_of_inlets, inlets)
    # @test length(inlets) == number_of_inlets; println("checkInletFiles - [OK]")
    #
    # @test_nowarn openBF.copyInputFilesToResultsFolder(project_name, inlets); println("copyInputFilesToResultsFolder - [1/5]")
    # @test isfile("$project_name.csv"); println("copyInputFilesToResultsFolder - [2/5]")
    # @test isfile("$project_name\_constants.yml"); println("copyInputFilesToResultsFolder - [3/5]")
    # @test isfile("$project_name\_inlet.dat"); println("copyInputFilesToResultsFolder - [4/5]")
    # @test isfile("$project_name\_2_inlet.dat"); println("copyInputFilesToResultsFolder - [OK]")
    #
    # model, model_header = openBF.readModelData("$project_name.csv")
    # @test typeof(model) == Array{Any,2}; println("readModelData - [1/2]")
    # @test model_header == ["name" "sn" "tn" "in" "L" "M" "Rp" "Rd" "E" "Pext" "R1" "R2" "C"]; println("readModelData - [OK]")
    #
    # constants = openBF.loadConstants("$project_name\_constants.yml")
    # @test typeof(constants) == Dict{Any,Any}; println("loadConstants - [OK]")
    #
    # @test openBF.checkConstants(constants) == constants; println("checkConstants - [OK]")
    #
    # blood = openBF.buildBlood(constants)
    # @test typeof(blood) == Blood; println("buildBlood - [OK]")
    #
    # inlet_data = openBF.loadInletData("$project_name\_inlet.dat")
    # @test typeof(inlet_data) == Array{Float64,2}; println("loadInletData - [OK]")
    #
    # heart = openBF.buildHeart(constants, inlet_data, 1)
    # @test typeof(heart) == Heart; println("buildHeart - [OK]")
    #
    # @inferred openBF.buildHearts(constants, inlets); println("buildHearts - [OK]")
    #
    # dx, invDx, halfDx = openBF.meshVessel(1.0, 10)
    # @test dx == 0.001; println("meshVessel - [1/3]")
    # @test invDx == 1000.0; println("meshVessel - [2/3]")
    # @test halfDx == 0.0005; println("meshVessel - [OK]")
    #
    # A0 = ones(Float64, 2)
    # gamma = copy(A0)
    # BCout = ["", "", "", ""]
    # BCout, outlet = openBF.detectCapillaries(BCout, blood, A0, gamma)
    # @test outlet == "none"; println("detectCapillaries - [1/20]")
    # @test BCout[1] == 0.0; println("detectCapillaries - [2/20]")
    # @test BCout[2] == 0.0; println("detectCapillaries - [3/20]")
    # @test BCout[3] == 0.0; println("detectCapillaries - [4/20]")
    # @test BCout[4] == 0.0; println("detectCapillaries - [5/20]")
    #
    # BCout = [1.0, "", "", ""]
    # BCout, outlet = openBF.detectCapillaries(BCout, blood, A0, gamma)
    # @test outlet == "reflection"; println("detectCapillaries - [6/20]")
    # @test BCout[1] == 1.0; println("detectCapillaries - [7/20]")
    # @test BCout[2] == 0.0; println("detectCapillaries - [8/20]")
    # @test BCout[3] == 0.0; println("detectCapillaries - [9/20]")
    # @test BCout[4] == 0.0; println("detectCapillaries - [10/20]")
    #
    # BCout = ["", "", 1.0, 1.0]
    # BCout, outlet = openBF.detectCapillaries(BCout, blood, A0, gamma)
    # @test outlet == "wk3"; println("detectCapillaries - [11/20]")
    # @test BCout[1] == 0.0; println("detectCapillaries - [12/20]")
    # @test isapprox(BCout[2],  1298, atol=1.0); println("detectCapillaries - [13/20]")
    # @test isapprox(BCout[3], -1297, atol=1.0); println("detectCapillaries - [14/20]")
    # @test BCout[4] == 1.0; println("detectCapillaries - [15/20]")
    #
    # BCout = ["", 1.0, 1.0, 1.0]
    # BCout, outlet = openBF.detectCapillaries(BCout, blood, A0, gamma)
    # @test outlet == "wk3"; println("detectCapillaries - [16/20]")
    # @test BCout[1] == 0.0; println("detectCapillaries - [17/20]")
    # @test BCout[2] == 1.0; println("detectCapillaries - [18/20]")
    # @test BCout[3] == 1.0; println("detectCapillaries - [19/20]")
    # @test BCout[3] == 1.0; println("detectCapillaries - [OK]")
    #
    # m_row = model[1,:]
    # @test_nowarn openBF.parseModelRow(m_row); println("parseModelRow - [OK]")
    # vessel_name, sn, tn, rn, L, M, Rp, Rd, E, Pext, BCout =  openBF.parseModelRow(m_row)
    #
    # vessel = openBF.buildVessel(1, m_row, heart, blood)
    # @test typeof(vessel) == Vessel; println("buildVessel - [OK]")
    #
    # vessels, edge_list = openBF.buildArterialNetwork(model, heart, blood)
    # @test typeof(vessels) == Array{Vessel,1}; println("buildArterialNetwork - [1/2]")
    # @test typeof(edge_list) == Array{Int,2}; println("buildArterialNetwork - [OK]")
    #
    # cd("..")
    # rm("$project_name\_results", recursive=true)
    #
    # inputs = openBF.loadSimulationFiles(project_name)
    # @test typeof(inputs[1]) == Dict{Any,Any}; println("loadSimulationFiles - [1/5]")
    # @test typeof(inputs[2]) == Array{Any,2}; println("loadSimulationFiles - [2/5]")
    # @test typeof(inputs[3]) == Array{Any,1}; println("loadSimulationFiles - [3/5]")
    # @test typeof(inputs[4]) == Blood; println("loadSimulationFiles - [4/5]")
    # @test typeof(inputs[5]) == Float64; println("loadSimulationFiles - [OK]")
    #
    # cd("..")
    # rm("$project_name\_results", recursive=true)

    @test_nowarn openBF.loadSimulationFiles("test.yml")

    data = openBF.loadYAMLFile("test.yml")
    @test_throws ErrorException openBF.loadYAMLFile("notest.yml")
    @test typeof(data) == Dict{Any,Any}

    @test_nowarn openBF.checkSections(data)
    @test_nowarn openBF.checkSection(data, "blood", ["mu", "rho"])
    @test_nowarn openBF.checkSection(data, "solver", ["Ccfl", "cycles", "jump", "convergence tollerance"])

    delete!(data, "project name")
    @test_throws ErrorException openBF.checkSections(data)

    delete!(data["blood"], "mu")
    @test_throws ErrorException openBF.checkSection(data, "blood", ["mu", "rho"])

    @test_nowarn openBF.checkNetwork(data["network"])
    delete!(data["network"][1], "inlet")
    @test_throws ErrorException openBF.checkNetwork(data["network"])

    delete!(data["network"][2], "R0")
    @test_throws ErrorException openBF.checkVessel(2, data["network"][2])

    data = openBF.loadYAMLFile("test.yml")
    blood = openBF.buildBlood(data["blood"])
    @test typeof(blood) == Blood
    @test blood.rho_inv == 1.0/blood.rho

    Rp, Rd = openBF.computeRadii(data["network"][1])
    @test Rp == 0.76e-2
    @test Rd == 0.7581e-2
    Rp, Rd = openBF.computeRadii(data["network"][2])
    @test Rp == Rd

    Pext = openBF.computePext(data["network"][1])
    @test Pext == 10000.0
    Pext = openBF.computePext(data["network"][2])
    @test Pext == 0.0

    M, dx, invDx, halfDx = openBF.meshVessel(data["network"][1], data["network"][1]["L"])
    @test M == data["network"][1]["M"]
    @test dx == 0.001
    @test isapprox(invDx, 1.0/0.001, atol=1e-3)
    @test halfDx == 0.0005
    M, dx, invDx, halfDx = openBF.meshVessel(data["network"][2], data["network"][2]["L"])
    @test M == 86
    M, dx, invDx, halfDx = openBF.meshVessel(data["network"][3], data["network"][3]["L"])
    @test M == 1000

    outlet, Rt, R1, R2, Cc = openBF.addOutlet(data["network"][1])
    @test outlet == "wk2"
    @test Rt == 0.0
    @test R1 == 0.0
    @test R2 == 6.8123e7
    @test Cc == 3.6664e-10
    outlet, Rt, R1, R2, Cc = openBF.addOutlet(data["network"][2])
    @test outlet == "reflection"
    @test Rt == 0.5
    @test R1 == 0.0
    @test R2 == 0.0
    @test Cc == 0.0
    outlet, Rt, R1, R2, Cc = openBF.addOutlet(data["network"][3])
    @test outlet == "wk3"
    @test Rt == 0.0
    @test R1 == 6.8123e7
    @test R2 == 3.1013e9
    @test Cc == 3.6664e-10

    @test isapprox(openBF.computeViscousTerm(data["network"][1], blood), 0.27, atol=1e-2)
    @test isapprox(openBF.computeViscousTerm(data["network"][2], blood), 0.27, atol=1e-2)

    inlet_data = openBF.loadInletData(data["network"][1]["inlet file"])
    @test inlet_data[1,2] != 0.0

    inlet, heart = openBF.buildHeart(data["network"][1])
    @test inlet
    @test heart.input_data[1,2] != 0.0
    @test heart.inlet_number == 1
    @test heart.inlet_type == "Q"
    inlet, heart = openBF.buildHeart(data["network"][2])
    @test ~inlet
    @test heart.input_data[1,2] == 0.0
    @test heart.inlet_number == 0
    @test heart.inlet_type == "none"

    vessel = openBF.buildVessel(1, data["network"][1], blood)
    @test typeof(vessel) == Vessel

    vessels, edges = openBF.buildArterialNetwork(data["network"], blood)
    @test length(vessels) == 3
    @test length(edges) == 12

    openBF.makeResultsFolder(data)
    cd("..")
    @test isdir("test_results")
end
