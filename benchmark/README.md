# benchmark

In the REPL


```julia
using PkgBenchmark
using openBF

bmresults = benchmarkpkg(openBF)
export_markdown("bm_results.md", bmresults)

jresults = judge(openBF, "v2.2.0") # or any other tag following this release
export_markdown("judge_bm_results.md", jresults)
```
