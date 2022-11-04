using BenchmarkTools
using openBF

const SUITE = BenchmarkGroup()

function benchmarkNetwork(network::String)
    cd("test/$(network)/")
    openBF.runSimulation("$(network).yml", verbose=false, out_files=true)
    rm("$(network)_results/", recursive=true)
    cd("../../")
end

for network in ("single-artery", "bifurcation", "tapering", "conjunction", "aspirator")
    SUITE[network] = @benchmarkable benchmarkNetwork($network)
end
