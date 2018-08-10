test_folder = "aspirator"
cd(test_folder)

openBF.runSimulation("$test_folder.yml", verbose=true, out_files=false)

p = readdlm("$test_folder"*"_results/d_P.last")

@test isapprox(minimum(p[:,end]), -1093.0, atol=1)
@test isapprox(mean(p[:,end]), 549.0, atol=1)
@test isapprox(maximum(p[:,end]), 5631.0, atol=1)

rm("$test_folder"*"_results/", recursive=true)
cd("..")
