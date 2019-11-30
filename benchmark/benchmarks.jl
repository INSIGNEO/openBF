using BenchmarkTools
using openBF

const SUITE = BenchmarkGroup()

A = pi*(0.5e-2)^2
A0 = pi*(0.45e-2)^2
sAoverA0 = sqrt(A/A0)
β = 1e6
Pext = 0.0

SUITE["pressure"] = BenchmarkGroup()
SUITE["pressure"]["normal"] = @benchmarkable openBF.pressure(A, A0, β, Pext)
SUITE["pressure"]["sqrt(A/A0)"] = @benchmarkable openBF.pressure(sAoverA0, β, Pext)

sA = sqrt(A)
γ = 1.0

SUITE["waveSpeed"] = BenchmarkGroup()
SUITE["waveSpeed"]["normal"] = @benchmarkable openBF.waveSpeed(A, γ)
SUITE["waveSpeed"]["sqrt(A)"] = @benchmarkable openBF.waveSpeedSA(sA, γ)