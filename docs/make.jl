using ToyClimaAtmos
using Documenter

DocMeta.setdocmeta!(ToyClimaAtmos, :DocTestSetup, :(using ToyClimaAtmos); recursive = true)

makedocs(;
    modules = [ToyClimaAtmos],
    authors = "Gabriele Bozzola <sbozzolator@gmail.com> and contributors",
    repo = "https://github.com/sbozzolo/ToyClimaAtmos.jl/blob/{commit}{path}#{line}",
    sitename = "ToyClimaAtmos.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)
