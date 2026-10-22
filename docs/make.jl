import Pkg
local_build = "--local" in ARGS
local_build || Pkg.update()
using RationalFunctionApproximation
using Documenter
using DocumenterCitations
local_build || @eval using DocumenterVitepress
bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

DocMeta.setdocmeta!(RationalFunctionApproximation, :DocTestSetup, :(using RationalFunctionApproximation); recursive=true)

makedocs(;
    modules=[RationalFunctionApproximation],
    authors="Toby Driscoll <driscoll@udel.edu>",
    repo=Remotes.GitHub("complexvariables", "RationalFunctionApproximation.jl"),
    sitename="RationalFunctionApproximation.jl",
    doctest=false,
    plugins=[bib],
    format=local_build ? Documenter.HTML() : DocumenterVitepress.MarkdownVitepress(;
        repo = "https://github.com/complexvariables/RationalFunctionApproximation.jl",
        devbranch = "main",
    ),
    pages=[
        "Introduction" => "index.md",
        "Installation" => "install.md",
        "Algorithms" => "algorithms.md",
        "Quadratic Thiele" => "qtcf.md",
        "Domains" => "domains.md",
        "Discrete data" => "discrete.md",
        "Minimax" => "minimax.md",
        "Usage from Python" => "python.md",
        "Functions" => "functions.md",
    ],
)

if !local_build
    DocumenterVitepress.deploydocs(;
        repo = "github.com/complexvariables/RationalFunctionApproximation.jl",
        target = joinpath(@__DIR__, "build"),
        branch = "gh-pages",
        devbranch = "main",
        push_preview = true,
    )
end
