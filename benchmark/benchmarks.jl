using openBF
using BenchmarkTools
using JSON

# ---------------------------------------------------------------------------
# PkgBenchmark-compatible SUITE
# Used by: benchmarkpkg(openBF) / judge(openBF, "vtag")
# ---------------------------------------------------------------------------

function _bmodel(dir::String, modelname::String)
    cwd = pwd()
    cd(dir)
    openBF.run_simulation("$modelname.yaml", verbose=false)
    cd(cwd)
end

const _ROOT = joinpath(@__DIR__, "..", "models")
const _B15  = joinpath(_ROOT, "boileau2015")
const _A07  = joinpath(_ROOT, "alastruey2007")

const SUITE = BenchmarkGroup()

SUITE["single_vessel"] = BenchmarkGroup()
SUITE["single_vessel"]["cca"] = @benchmarkable _bmodel(joinpath($_B15, "cca"), "cca")
SUITE["single_vessel"]["uta"] = @benchmarkable _bmodel(joinpath($_B15, "uta"), "uta")

SUITE["bifurcation"] = BenchmarkGroup()
SUITE["bifurcation"]["ibif"] = @benchmarkable _bmodel(joinpath($_B15, "ibif"), "ibif")

SUITE["circulation"] = BenchmarkGroup()
SUITE["circulation"]["adan56"]           = @benchmarkable _bmodel(joinpath($_B15, "adan56"), "adan56")
SUITE["circulation"]["circle_of_willis"] = @benchmarkable _bmodel($_A07, "circle_of_willis")

# ---------------------------------------------------------------------------
# Manual regression harness
# Usage:
#   include("benchmark/benchmarks.jl")
#   run_suite("after_my_change")
#   judge_against("after_my_change", "baseline")
# ---------------------------------------------------------------------------

const HARNESS_MODELS = [
    ("cca",              joinpath(_B15, "cca")),
    ("ibif",             joinpath(_B15, "ibif")),
    ("adan56",           joinpath(_B15, "adan56")),
    ("circle_of_willis", _A07),
]

const RESULTS_PATH = joinpath(@__DIR__, "results")
isdir(RESULTS_PATH) || mkpath(RESULTS_PATH)

function bench_model(name::String, dir::String)
    cwd = pwd()
    cd(dir)
    t = @benchmark openBF.run_simulation($name * ".yaml", verbose=false) seconds=60 samples=5
    cd(cwd)
    t
end

function run_suite(label::String)
    results = Dict{String,Any}()
    for (m, sub) in HARNESS_MODELS
        println("Benchmarking $m ...")
        b = bench_model(m, sub)
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
    for (m, _) in HARNESS_MODELS
        (haskey(a, m) && haskey(b, m)) || continue
        Δ = (a[m]["min_time_ns"] - b[m]["min_time_ns"]) / b[m]["min_time_ns"] * 100
        @info "$m" time_change_pct=round(Δ, digits=1) allocs_before=b[m]["allocs"] allocs_after=a[m]["allocs"]
    end
end
