push!(LOAD_PATH, "../src/")
using Documenter, openBF

makedocs(
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename = "openBF",
    pages = [
        "Home" => "index.md",
        "Overview" => "man/overview.md",
        "Quickstart" => "man/quickstart.md",
        "Configuration" => "man/config.md",
        "Examples" => "man/examples.md",
        "Plotting" => "man/plotting.md",
        "Benchmarking" => "man/benchmarking.md",
    ],
)

deploydocs(
    repo = "github.com/INSIGNEO/openBF.git",
)
