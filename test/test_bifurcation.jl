test_folder = "bifurcation"
cd(test_folder)

openBF.runSimulation("$test_folder.yml", verbose=true, out_files=false)

p = readdlm("$test_folder"*"_results/d1_P.last")

@test isapprox(minimum(p[:,end]), 8827.0, atol=20)
@test isapprox(mean(p[:,end]), 12622.0, atol=20)
@test isapprox(maximum(p[:,end]), 17695.0, atol=20)

rm("$test_folder"*"_results/", recursive=true)
cd("..")

