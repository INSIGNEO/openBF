using Test
using DelimitedFiles
using Statistics
using openBF

@testset "openBF.jl" begin

    #unit tests
    # println("Test initialise.jl functions")
    # include("test_initialise.jl")

    # println("Test boundary_conditions.jl functions")
    # include("test_boundary_conditions.jl")

    # println("Test solver.jl functions")
    # include("test_solver.jl")

    #integration tests
    # include("test_single-artery.jl")
    include("test_conjunction.jl")
    include("test_bifurcation.jl")
    include("test_aspirator.jl")
    include("test_tapering.jl")
end
