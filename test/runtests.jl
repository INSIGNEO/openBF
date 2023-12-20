using Test
using DelimitedFiles
using Statistics
using openBF

@testset "openBF.jl" begin
    include("test_network.jl")
    include("test_adan.jl")
end
