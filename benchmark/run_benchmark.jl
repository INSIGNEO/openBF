using PkgBenchmark
bmr = benchmarkpkg("openBF")
export_markdown("bm_results.md", bmr)