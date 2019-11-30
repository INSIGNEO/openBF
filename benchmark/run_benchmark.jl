using PkgBenchmark
bmr = benchmarkpkg("openBF")
export_markdown("bm_results.md", bmr)

jbm = judge("openBF", "9cfab60ff269e43a5ed752c6ef3d6fdaf132efa1")
export_markdown("bm_comparison.md", jbm)
