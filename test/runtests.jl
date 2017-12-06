using Base.Test

function runTest(test_folder)

  cd(test_folder)
  try
    run(`cp ../../main.jl main.jl`)
    # run(`cp ../../src/BTypes.jl BTypes.jl`)
    run(`julia main.jl $test_folder`)
    run(`rm BTypes.jl`)
    cd("..")
    return true

  catch
    # run(`rm BTypes.jl`)
    cd("..")
    return false
  end
end

# Run tests

tic()
println("\nTest 1 - Single artery")
@time @test runTest("single-artery")

println("\nTest 2 - Bifurcation")
@time @test runTest("bifurcation")

println("\nTest 3 - External pressure")
@time @test runTest("external-pressure")
toc()
