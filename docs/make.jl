push!(LOAD_PATH, ".../src/")
using Documenter

makedocs(sitename="caitlyn.jl",
        pages = [
            "Home" => "index.md",
            "Algorithm" => "algorithm.md",
            "Usage" => "usage.md"
        ])
