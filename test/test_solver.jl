data = openBF.loadSimulationFiles("test.yml")
openBF.makeResultsFolder(data)
blood = openBF.buildBlood(data["blood"])
vessels, edges = openBF.buildArterialNetwork(data["network"], blood, 100)

Ccfl = data["solver"]["Ccfl"]
heart = vessels[1].heart

A = vessels[1].A[1]
A0 = vessels[1].A0[1]
beta = vessels[1].beta[1]
gamma = vessels[1].gamma[1]
Pext = vessels[1].Pext

p = openBF.pressure(A, A0, beta, Pext)
@test p == 1e4
@test p == openBF.pressure(sqrt(A/A0), beta, Pext)

c = openBF.waveSpeed(A, gamma)
@test isapprox(c, 6.33, atol=1e-2)

dt = openBF.calculateDeltaT(vessels, Ccfl)
@test isapprox(dt, 1e-4, atol=1e-4)

@test_nowarn openBF.solveModel(vessels, edges, blood, dt, 0.0)

q = vessels[1].Q[1]
p = vessels[1].P[end]
@test_nowarn openBF.solveVessel(vessels[1], blood, dt, 0.0)
@test q != vessels[1].Q[1]
@test p != vessels[1].P[end]

@test_nowarn openBF.solveOutlet(1, vessels[1], blood, vessels, edges, dt)

v = vessels[1]
@test_nowarn openBF.computeLimiter(v, v.vA, v.invDx, v.dU, v.slopesA)
@test_nowarn openBF.computeFlux(v, v.Al, v.Ql, v.Fl)

@test openBF.maxMod(1.0, 0.0) == 1.0
@test openBF.maxMod(0.0, 0.5) == 0.5

@test openBF.minMod(-1.0, 1.0) == 0.0
@test openBF.minMod(1.0, -1.0) == 0.0
@test openBF.minMod(0.5, 1.0) == 0.5
@test openBF.minMod(1.0, 0.5) == 0.5

@test_nowarn openBF.superBee(v, v.dU, v.slopesA)

cd("..")
rm("test_results", recursive=true)
