using openBF
using YAML
using DelimitedFiles
using Serialization
using Statistics
using Printf

const MODEL_PATHS = Dict(
    "cca"              => joinpath(@__DIR__, "..", "models", "boileau2015", "cca"),
    "ibif"             => joinpath(@__DIR__, "..", "models", "boileau2015", "ibif"),
    "adan56"           => joinpath(@__DIR__, "..", "models", "boileau2015", "adan56"),
    "circle_of_willis" => joinpath(@__DIR__, "..", "models", "alastruey2007"),
)

"""
    capture_waveforms(model_name) -> (data, config)

Run a simulation and return the last-cycle waveform matrices keyed by vessel label.
`data[label][field]` is a Matrix with columns [t, val_node1, val_node2, val_node3, val_node4, val_nodeM].
"""
function capture_waveforms(model_name::String)
    dir = abspath(MODEL_PATHS[model_name])
    yaml = joinpath(dir, "$model_name.yaml")
    savedir = mktempdir(; prefix="openbf_$(model_name)_", cleanup=false)

    cwd = pwd()
    openBF.run_simulation(yaml, verbose=false, savedir=savedir)
    cd(cwd)

    config = YAML.load_file(yaml)
    data = Dict{String, Dict{String, Matrix{Float64}}}()
    for v in config["network"]
        label = v["label"]
        waves = Dict{String, Matrix{Float64}}()
        for field in ("Q", "P", "A", "u")
            fpath = joinpath(savedir, "$(label)_$field.last")
            isfile(fpath) && (waves[field] = readdlm(fpath, Float64))
        end
        !isempty(waves) && (data[label] = waves)
    end

    rm(savedir; recursive=true, force=true)
    data, config
end

"""
    write_reference(model_name)

Capture and serialize waveforms to `test/ref_<model_name>.jls`.
"""
function write_reference(model_name::String)
    data, _ = capture_waveforms(model_name)
    path = joinpath(@__DIR__, "ref_$(model_name).jls")
    serialize(path, data)
    println("Written: $path")
    data
end

"""
    load_reference(model_name) -> data dict
"""
load_reference(model_name::String) =
    deserialize(joinpath(@__DIR__, "ref_$(model_name).jls"))

const ALL_MODELS = ("cca", "ibif", "adan56", "circle_of_willis")

"""
    check_all(; rtol=0, atol=0) -> Bool

Run all four models, compare waveforms against stored references.
Returns true if every vessel in every model is within tolerance.
"""
function check_all(; rtol=0.0, atol=0.0)
    all_ok = true
    for m in ALL_MODELS
        ref = load_reference(m)
        data, _ = capture_waveforms(m)
        diffs = waveform_diff(ref, data; field="P")
        max_diff = isempty(diffs) ? 0.0 : maximum(values(diffs))
        tol = rtol + atol
        ok = max_diff <= tol
        all_ok &= ok
        status = ok ? "PASS" : "FAIL"
        @printf "  [%s] %s  max_rel_diff=%.6f  (tol=%.6f)\n" status m max_diff tol
    end
    all_ok
end

"""
    waveform_diff(data_a, data_b; field="P") -> Dict{label => max_rel_diff}

Per-vessel maximum pointwise relative difference on `field` (columns 2:end, skipping time).
"""
function waveform_diff(data_a, data_b; field="P")
    result = Dict{String, Float64}()
    for label in intersect(keys(data_a), keys(data_b))
        wa = get(data_a[label], field, nothing)
        wb = get(data_b[label], field, nothing)
        (wa === nothing || wb === nothing) && continue
        a = wa[:, 2:end]
        b = wb[:, 2:end]
        denom = maximum(abs, a)
        result[label] = denom > 0 ? maximum(abs, a .- b) / denom : 0.0
    end
    result
end
