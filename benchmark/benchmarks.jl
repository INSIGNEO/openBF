using BenchmarkTools
using openBF

cd("test/")
data = openBF.loadSimulationFiles("test.yml")
blood = openBF.buildBlood(data["blood"])
vessels, edges = openBF.buildArterialNetwork(data["network"], blood, 100)

const SUITE = BenchmarkGroup()
SUITE["solver"] = BenchmarkGroup()

A = vessels[1].A[1]
A0 = vessels[1].A0[1]
β = vessels[1].beta[1]
Pext = vessels[1].Pext
SUITE["solver"]["pressure"] = @benchmarkable openBF.pressure(A, A0, β, Pext)

sAoverA0 = sqrt(A/A0)
SUITE["solver"]["pressure_opt"] = @benchmarkable openBF.pressure(sAoverA0, β,
	Pext)

γ = vessels[1].gamma[1]
SUITE["solver"]["waveSpeed"] = @benchmarkable openBF.waveSpeed(A, γ)

sA = sqrt(A)
SUITE["solver"]["waveSpeedSA"] = @benchmarkable openBF.waveSpeedSA(sA, γ)

Ccfl = data["solver"]["Ccfl"]
SUITE["solver"]["calculateDeltaT"] = @benchmarkable openBF.calculateDeltaT(vessels,
	Ccfl)
dt = openBF.calculateDeltaT(vessels, Ccfl)


SUITE["solver"]["solveModel"] = @benchmarkable openBF.solveModel(vessels, edges,
	blood, dt, 0.0)

SUITE["solver"]["solveVessel"] = @benchmarkable openBF.solveVessel(vessels[1],
	blood, dt, 0.0)

# SUITE["solver"]["solveOutlet"] = @benchmarkable openBF.solveOutlet(1,
# 	vessels[1], blood, vessels, edges, dt)

v = vessels[1]
SUITE["solver"]["computeLimiter"] = @benchmarkable openBF.computeLimiter(v,
	v.vA, v.invDx, v.dU, v.slopesA)
SUITE["solver"]["computeFlux"] = @benchmarkable openBF.computeFlux(v, v.Al,
	v.Ql, v.Fl)

SUITE["solver"]["maxMod"] = @benchmarkable openBF.maxMod(1.0, 0.0)
SUITE["solver"]["minMod"] = @benchmarkable openBF.minMod(-1.0, 1.0) == 0.0

SUITE["solver"]["superBee"] = @benchmarkable openBF.superBee(v, v.dU, v.slopesA)
