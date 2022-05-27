using Documenter
using PeriodicGraphEmbeddings, PeriodicGraphs, Graphs

DocMeta.setdocmeta!(PeriodicGraphEmbeddings, :DocTestSetup, quote
    using PeriodicGraphEmbeddings, PeriodicGraphs, Graphs
end; recursive=true)

makedocs(
    sitename = "PeriodicGraphEmbeddings.jl",
    format = Documenter.HTML(),
    modules = [PeriodicGraphEmbeddings],
    pages = [
        "Home" => "index.md",
        "Types" => "types.md",
        "Symmetries" => "symmetries.md",
        "I/O" => "io.md",
        "Utilities"  => "utilities.md"
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/Liozou/PeriodicGraphEmbeddings.jl.git"
)
