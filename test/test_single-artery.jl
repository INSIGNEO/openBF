test_folder = "single-artery"
cd(test_folder)

openBF.runSimulation("$test_folder.yml", verbose=true, out_files=false)

p = readdlm("$test_folder"*"_results/A1_P.last")

@test isapprox(minimum(p[:,end]), 9522.0, atol=20)
@test isapprox(mean(p[:,end]), 12712.0, atol=10)
@test isapprox(maximum(p[:,end]), 16772.0, atol=10)

rm("$test_folder"*"_results/", recursive=true)
cd("..")
