using Base.Test

function runTest(test_folder)

  cd(test_folder)
  try
    run(`julia6 main.jl`)
    cd("..")
    return true

  catch
    cd("..")
    return false
  end
end

# Run tests

tic()
println("\nTest 1")
@time @test runTest("single-artery")

println("\nTest 2")
@time @test runTest("bifurcation")

println("\nTest 3")
@time @test runTest("external-pressure")
toc()
