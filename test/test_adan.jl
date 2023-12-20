println("test adan")
test_folder = "adan"
cd(test_folder)

rm("adan_results/", recursive = true, force = true)

@test_nowarn openBF.run_simulation("adan56.yaml", verbose = false, out_files = false)

for artery in [
    "aortic_arch_I",
    "abdominal_aorta_V",
    "common_carotid_R",
    "common_iliac_R",
    "radial_R",
    "posterior_interosseous_R",
]

    P = readdlm("adan56_results/$(artery)_P.last")
    p = readdlm("ref/$(artery)_P.last")
    Q = readdlm("adan56_results/$(artery)_Q.last")
    q = readdlm("ref/$(artery)_Q.last")

    @test isapprox(P[:,4], p[:,4], atol=1.0)
    @test isapprox(Q[:,4], q[:,4], atol=1e-6)
end

rm("adan56_results/", recursive = true)
cd("..")
