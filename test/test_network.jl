test_folder = "network"
cd(test_folder)

rm("$test_folder" * "_results/", recursive = true, force = true)

@test_nowarn openBF.run_simulation("$test_folder.yaml", verbose = false, out_files = false)

rm("$test_folder" * "_results/", recursive = true)
cd("..")
