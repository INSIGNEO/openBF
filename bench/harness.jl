using openBF
using BenchmarkTools
using JSON

const MODELS = ["cca", "ibif", "adan56"]
const RESULTS_PATH = joinpath(@__DIR__, "results")
isdir(RESULTS_PATH) || mkpath(RESULTS_PATH)

function bench_model(name::String)
    cwd = pwd()
    cd(joinpath(@__DIR__, "..", "models", "boileau2015", name))
    t = @benchmark openBF.run_simulation($name * ".yaml", verbose=false) seconds=60 samples=5
    cd(cwd)
    t
end

function run_suite(label::String)
    results = Dict{String,Any}()
    for m in MODELS
        println("Benchmarking $m ...")
        b = bench_model(m)
        results[m] = Dict(
            "min_time_ns"    => minimum(b).time,
            "median_time_ns" => median(b).time,
            "mem_bytes"      => minimum(b).memory,
            "allocs"         => minimum(b).allocs,
            "gc_pct"         => minimum(b).gctime / minimum(b).time,
        )
        println("  $m: min=$(round(minimum(b).time/1e9, digits=3))s  allocs=$(minimum(b).allocs)")
    end
    open(joinpath(RESULTS_PATH, "$label.json"), "w") do io
        JSON.print(io, results, 2)
    end
    results
end

function judge_against(label::String, baseline::String="baseline")
    a = JSON.parsefile(joinpath(RESULTS_PATH, "$label.json"))
    b = JSON.parsefile(joinpath(RESULTS_PATH, "$baseline.json"))
    for m in MODELS
        t_new = a[m]["min_time_ns"]
        t_old = b[m]["min_time_ns"]
        Δ = (t_new - t_old) / t_old * 100
        alloc_new = a[m]["allocs"]
        alloc_old = b[m]["allocs"]
        @info "$m" time_change_pct=round(Δ, digits=1) allocs_before=alloc_old allocs_after=alloc_new
    end
end
