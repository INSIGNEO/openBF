# benchmark

## PkgBenchmark workflow

```julia
using PkgBenchmark, openBF

# run full suite against current code
results = benchmarkpkg(openBF)
export_markdown("bm_results.md", results)

# compare against a tagged release
judge_results = judge(openBF, "v2.2.0")
export_markdown("judge_bm_results.md", judge_results)
```

## Manual regression harness

```julia
include("benchmark/benchmarks.jl")

# run and save results
run_suite("after_my_change")

# compare against a saved baseline
judge_against("after_my_change", "after_step4_26")
```

Results are saved as JSON files in `benchmark/results/`.
