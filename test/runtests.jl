using Base.Test
using openBF

function runTest(test_folder)
    println(" ")
    cd(test_folder)
    try
        openBF.runSimulation("$test_folder.yml", verbose=true, out_files=false)
        cd("..")
    return true

    catch
        cd("..")
        return false
    end
end

@testset "openBF.jl" begin
    println("Test initialise.jl functions")
    include("test_initialise.jl")

    println("Test boundary_conditions.jl functions")
    include("test_boundary_conditions.jl")

    println("Test solver.jl functions")
    include("test_solver.jl")

    @time @test runTest("single-artery")
    @time @test runTest("conjunction")
    @time @test runTest("bifurcation")
    @time @test runTest("aspirator")

    # println("\nTest External pressure")
    # @time @test runTest("external-pressure")
    #
    # println("\nTest Anastomosis")
    # @time @test runTest("anastomosis")
end
