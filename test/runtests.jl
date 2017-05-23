using Base.Test

function runTest(test_folder)

  cd(test_folder)
  try
    run(`julia main.jl`)
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
@time @test runTest("artery-vein")

println("\nTest 3")
@time @test runTest("bifurcation")
toc()
