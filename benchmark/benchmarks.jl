using openBF
using BenchmarkTools

function bmodel(modelname)
    cwd = pwd()
    println(cwd)
    cd(joinpath("models", "boileau2015", modelname))
    openBF.run_simulation("$modelname.yaml")
    cd(cwd)
end

const SUITE = BenchmarkGroup()

SUITE["single_vessel"] = BenchmarkGroup()
SUITE["single_vessel"]["cca"] = @benchmarkable bmodel("cca")
SUITE["single_vessel"]["uta"] = @benchmarkable bmodel("uta")

SUITE["bifurcation"] = BenchmarkGroup()
SUITE["bifurcation"]["ibif"] = @benchmarkable bmodel("ibif")
