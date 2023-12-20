using Test
using DelimitedFiles
using Statistics
using openBF

@testset "openBF.jl" begin
    println("Test network")
    include("test_network.jl")
    include("test_adan.jl")
end
