# Fast A/B: ibif has a bifurcation junction and well-configured windkessel outlets.
# Conjunction and anastomosis are covered at solver level by test_junction_k2/k3.
# Verifies use_generic_junctions=true matches legacy to 1e-10 relative.

@testset "Junction generic solver (A/B on ibif)" begin
    dir    = abspath(joinpath(@__DIR__, "..", "models", "boileau2015", "ibif"))
    yaml   = joinpath(dir, "ibif.yaml")
    labels = ["parent", "d1", "d2"]
    fields = ("Q", "P", "A", "u")
    tol    = 1e-10

    dir_leg = mktempdir(; prefix="openbf_leg_", cleanup=false)
    dir_gen = mktempdir(; prefix="openbf_gen_", cleanup=false)
    cwd     = pwd()

    openBF.run_simulation(yaml; verbose=false, savedir=dir_leg)
    cd(cwd)
    openBF.run_simulation(yaml; verbose=false, savedir=dir_gen, use_generic_junctions=true)
    cd(cwd)

    for label in labels, field in fields
        fleg = joinpath(dir_leg, "$(label)_$field.last")
        fgen = joinpath(dir_gen, "$(label)_$field.last")
        (isfile(fleg) && isfile(fgen)) || continue
        leg   = readdlm(fleg, Float64)[:, 2:end]
        gen   = readdlm(fgen, Float64)[:, 2:end]
        denom = maximum(abs, leg)
        denom > 0 && @test maximum(abs, leg .- gen) / denom ≤ tol
    end

    rm(dir_leg; recursive=true, force=true)
    rm(dir_gen; recursive=true, force=true)
end
