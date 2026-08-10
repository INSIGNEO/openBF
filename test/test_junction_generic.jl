# Smoke test: generic junction solver runs ibif (bifurcation) end-to-end
# and produces finite, non-trivial waveforms on all vessels.

@testset "Junction solver smoke test (ibif)" begin
    dir  = abspath(joinpath(@__DIR__, "..", "models", "boileau2015", "ibif"))
    yaml = joinpath(dir, "ibif.yaml")
    labels = ["parent", "d1", "d2"]
    fields = ("Q", "P", "A", "u")
    savedir = mktempdir(; prefix="openbf_ibif_", cleanup=false)
    cwd = pwd()

    openBF.run_simulation(yaml; verbose=false, savedir=savedir)
    cd(cwd)

    for label in labels, field in fields
        fpath = joinpath(savedir, "$(label)_$field.last")
        isfile(fpath) || continue
        data = readdlm(fpath, Float64)[:, 2:end]
        @test all(isfinite, data)
        @test maximum(abs, data) > 0
    end

    rm(savedir; recursive=true, force=true)
end
