data = openBF.loadSimulationFiles("test.yml")
# openBF.makeResultsFolder(data)
blood = openBF.buildBlood(data["blood"])
vessels, edges = openBF.buildArterialNetwork(data["network"], blood, 100)

Ccfl = data["solver"]["Ccfl"]
heart = vessels[1].heart

Pref = vessels[1].P[1]
Qref = vessels[1].Q[1]
dt = openBF.calculateDeltaT(vessels, Ccfl)

@test heart.inlet_type == "Q"
@test_nowarn openBF.setInletBC(0.0, dt, vessels[1])
@test vessels[1].Q[1] != Qref
@test vessels[1].P[1] == Pref

@test_nowarn openBF.inletCompatibility(dt, vessels[1], heart)
@test isapprox(vessels[1].A[1], 1.8e-4, atol=1e-5)
@test isapprox(vessels[1].P[1], 1e4, atol=1)

W1, W2 = openBF.riemannInvariants(1, vessels[1])
@test isapprox(W1, -W2, atol=1)

u, c = openBF.inverseRiemannInvariants(W1, W2)
@test isapprox(u, vessels[1].u[1], atol=1e-4)
@test c == vessels[1].c[1]

newA = openBF.areaFromPressure(vessels[1].P[1], vessels[1].A0[1], vessels[1].beta[1],
                                vessels[1].Pext)
@test newA == vessels[1].A[1]

@test isapprox(openBF.inputFromData(0.0, heart), -5.239e-07, atol=1e-7)

@test_nowarn openBF.setOutletBC(dt, vessels[3])
@test isapprox(vessels[3].A[end], 9.5e-5, atol=1e-5)
@test vessels[3].u[end] == 0.0

@test_nowarn openBF.threeElementWindkessel(dt, vessels[1])

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

@test_nowarn openBF.updateGhostCells(vessels)

# cd("..")
# rm("test_results", recursive=true)
