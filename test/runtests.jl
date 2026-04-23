using Test
using DelimitedFiles
using Statistics
using openBF

@testset "openBF.jl" begin
    include("test_network.jl")
    include("test_junction_type.jl")
    include("test_junction_k2.jl")
    include("test_junction_k3_bifurcation.jl")
    include("test_junction_k3_anastomosis.jl")
    include("test_golden_master.jl")
end
