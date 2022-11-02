@test_nowarn openBF.loadSimulationFiles("test.yml")

data = openBF.loadYAMLFile("test.yml")
@test_throws ErrorException openBF.loadYAMLFile("notest.yml")
@test typeof(data) == Dict{Any,Any}

@test_nowarn openBF.checkInputFile(data)

@test_nowarn openBF.checkSections(data)
@test_nowarn openBF.checkSection(data, "blood", ["mu", "rho"])
@test_nowarn openBF.checkSection(data, "solver", ["Ccfl", "cycles", "convergence tolerance"])

delete!(data, "project name")
@test_throws ErrorException openBF.checkSections(data)

delete!(data["blood"], "mu")
@test_throws ErrorException openBF.checkSection(data, "blood", ["mu", "rho"])

@test_nowarn openBF.checkNetwork(data["network"])

two_inlets_same_node = deepcopy(data)
two_inlets_same_node["network"][2]["inlet"] = "Q"
two_inlets_same_node["network"][2]["sn"] = 1
@test_throws ErrorException openBF.checkNetwork(two_inlets_same_node["network"])

delete!(data["network"][1], "inlet")
@test_throws ErrorException openBF.checkNetwork(data["network"])

delete!(data["network"][2], "R0")
@test_throws ErrorException openBF.checkVessel(2, data["network"][2])

data["network"][2]["sn"] = data["network"][2]["tn"]
@test_throws ErrorException openBF.checkVessel(2, data["network"][2])

data = openBF.loadYAMLFile("test_connectivity.yml")
@test_throws ErrorException openBF.checkNetwork(data["network"])

data["network"][3]["sn"] = 3
data["network"][4]["sn"] = 4
@test_throws ErrorException openBF.checkNetwork(data["network"])

data["network"][4]["outlet"] = "wk2"
data["network"][3]["sn"] = 6
@test_throws ErrorException openBF.checkNetwork(data["network"])

data = openBF.loadYAMLFile("test.yml")
blood = openBF.buildBlood(data["blood"])
@test typeof(blood) == Blood
@test blood.rho_inv == 1.0/blood.rho

Rp, Rd = openBF.computeRadii(data["network"][1])
@test Rp == 0.76e-2
@test Rd == 0.7581e-2
Rp, Rd = openBF.computeRadii(data["network"][2])
@test Rp == Rd

Pext = openBF.getPext(data["network"][1])
@test Pext == 10000.0
Pext = openBF.getPext(data["network"][2])
@test Pext == 0.0

M, dx, invDx, halfDx, invDxSq = openBF.meshVessel(data["network"][1], data["network"][1]["L"])
@test M == data["network"][1]["M"]
@test dx == 0.001
@test isapprox(invDx, 1.0/0.001, atol=1e-3)
@test isapprox(invDxSq, (1.0/0.001)^2, atol=1e-3)
@test halfDx == 0.0005
h0 = openBF.initialiseThickness(data["network"][1], M)
@test isequal(h0, 0.0)

M, dx, invDx, halfDx = openBF.meshVessel(data["network"][2], data["network"][2]["L"])
@test M == 86
h0 = openBF.initialiseThickness(data["network"][2], M)
@test ~isequal(h0, 0.0)

M, dx, invDx, halfDx, invDxSq = openBF.meshVessel(data["network"][3], data["network"][3]["L"])
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

@test isapprox(openBF.computeViscousTerm(data["network"][1], blood), 0.27*blood.rho_inv, atol=1e-2)
@test isapprox(openBF.computeViscousTerm(data["network"][2], blood), 0.27*blood.rho_inv, atol=1e-2)

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

vessel = openBF.buildVessel(1, data["network"][1], blood, 100)
@test typeof(vessel) == Vessel

vessels, edges = openBF.buildArterialNetwork(data["network"], blood, 100)
@test length(vessels) == 3
@test length(edges) == 9

v = vessels[1]
R1, R2 = openBF.computeWindkesselInletImpedance(v.R2, blood, v.A0, v.gamma)
@test R1 + R2 == v.R2

openBF.makeResultsFolder(data, "test.yml")
cd("..")
@test isdir("test_results")
rm("test_results", recursive=true)
