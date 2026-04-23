# Benchmarking

openBF ships a benchmark suite in `benchmark/` compatible with [PkgBenchmark.jl](https://github.com/JuliaCI/PkgBenchmark.jl) and a lightweight manual harness for tracking regressions across changes.

## PkgBenchmark workflow

```julia
using PkgBenchmark, openBF

# run the full suite against the current code
results = benchmarkpkg(openBF)
export_markdown("bm_results.md", results)

# compare against a tagged release
judge_results = judge(openBF, "v2.2.0")
export_markdown("judge_bm_results.md", judge_results)
```

The suite covers five models grouped by complexity:

| Group | Models |
|-------|--------|
| `single_vessel` | `cca`, `uta` |
| `bifurcation` | `ibif` |
| `circulation` | `adan56`, `circle_of_willis` |

## Manual regression harness

For tracking performance across commits without PkgBenchmark:

```julia
include("benchmark/benchmarks.jl")

# run all models and save results to benchmark/results/<label>.json
run_suite("after_my_change")

# compare two saved runs
judge_against("after_my_change", "baseline")
```

`judge_against` prints the percentage change in minimum time and allocation count for each model. Results are stored locally in `benchmark/results/` (not tracked by git).

## Benchmarked models

The harness runs `cca`, `ibif`, `adan56`, and `circle_of_willis`. Each model is sampled 5 times with a 60-second budget; the minimum time is reported. Allocation counts are a more stable signal than wall-clock time on shared machines.
