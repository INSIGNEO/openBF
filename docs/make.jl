using Documenter, openBF

makedocs(
    format = :html,
    sitename = "openBF",
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
            "Guide" => "man/guide.md",
            "man/config.md",
            "man/examples.md"
        ]
    ])
