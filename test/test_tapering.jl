test_folder = "tapering"
cd(test_folder)

openBF.runSimulation("$test_folder.yml", verbose=true, out_files=false)

p = readdlm("$test_folder"*"_results/A1_u.last")

@test isapprox(sum(p[:,2]), 14.65, atol=1)
@test isapprox(sum(p[:,4]), 18.0, atol=1)
@test isapprox(sum(p[:,end]), 23.4, atol=1)

rm("$test_folder"*"_results/", recursive=true)
cd("..")
