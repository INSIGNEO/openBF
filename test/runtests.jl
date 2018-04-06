using Base.Test
using openBF

function runTest(test_folder)
    println(" ")
    cd(test_folder)
    try
        openBF.runSimulation("$test_folder.yml", verbose=true, out_files=false)
        rm("$test_folder\_results", recursive=true)
        cd("..")
    return true

    catch
        cd("..")
        return false
    end
end

@testset "openBF.jl" begin
    #unit tests
    println("Test initialise.jl functions")
    include("test_initialise.jl")

    println("Test boundary_conditions.jl functions")
    include("test_boundary_conditions.jl")

    println("Test solver.jl functions")
    include("test_solver.jl")

    #integration tests
    @test runTest("single-artery")
    @test runTest("conjunction")
    @test runTest("bifurcation")
    @test runTest("aspirator")
end
