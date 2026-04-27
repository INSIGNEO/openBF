#!/usr/bin/env julia
# Usage:
#   julia run.jl <yaml>              # run simulation
#   julia run.jl --validate <yaml>  # validate network topology only
using openBF

if "--validate" in ARGS
    yaml = ARGS[findfirst(a -> !startswith(a, "-"), ARGS)]
    ok   = openBF.validate_network(yaml)
    exit(ok ? 0 : 1)
elseif !isempty(ARGS)
    openBF.run_simulation(ARGS[1])
else
    println("Usage: julia run.jl [--validate] <network.yaml>")
    exit(1)
end
