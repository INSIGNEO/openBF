using Base.Test
using openBF

function runTest(test_folder)

  cd(test_folder)
  try
    run(`cp ../../main.jl main.jl`)
    run(`$JULIA_HOME/julia main.jl $test_folder`)
    cd("..")
    return true

  catch
    cd("..")
    return false
  end
end

# Run tests

# println("\nTest 1 - Single artery")
# @time @test runTest("single-artery")
#
# println("\nTest 2 - Bifurcation")
# @time @test runTest("bifurcation")
#
# println("\nTest 3 - External pressure")
# @time @test runTest("external-pressure")
#
# println("\nTest 4 - Anastomosis")
# @time @test runTest("anastomosis")
#
# println("\nTest 5 - Conjunction")
# @time @test runTest("conjunction")
#
# println("\nTest 5 - Aspirator")
# @time @test runTest("aspirator")
@testset "openBF.jl" begin
    include("test_initialise.jl")
    include("test_boundary_conditions.jl")
    include("test_solver.jl")
end
