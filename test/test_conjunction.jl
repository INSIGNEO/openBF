test_folder = "conjunction"
cd(test_folder)

openBF.runSimulation("$test_folder.yml", verbose=true, out_files=false)

p = readdlm("$test_folder"*"_results/d1_P.last")

@test isapprox(minimum(p[:,end]), -1771.0, atol=1)
@test isapprox(mean(p[:,end]), 623.0, atol=1)
@test isapprox(maximum(p[:,end]), 6529.0, atol=1)

rm("$test_folder"*"_results/", recursive=true)
cd("..")
