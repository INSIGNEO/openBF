using BenchmarkTools
using openBF

const SUITE = BenchmarkGroup()

cd("test/single-artery")
SUITE["single_artery"] = BenchmarkGroup()
SUITE["single_artery"]["run"] = @benchmarkable openBF.runSimulation("single-artery.yml", out_files=false)
rm("single-artery_results/", recursive=true)

