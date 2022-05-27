"""
PeriodicGraphEmbeddings

A module for manipulating embeddings of periodic graphs in euclidean space.

Symmetry properties are provided on graphs up to dimension 3 through
[spglib](https://spglib.github.io/spglib/).
"""
module PeriodicGraphEmbeddings

using PeriodicGraphs, StaticArrays, Graphs
import LinearAlgebra: det, norm, istriu, cross, dot

import Base: ==
using Tokenize

include("utils.jl")
include("types.jl")
include("other.jl")
include("symmetries.jl")
include("io.jl")

include("precompile.jl")

end
