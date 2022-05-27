# Types

## Manual

PeriodicGraphEmbeddings provide the new type [`PeriodicGraphEmbedding`](@ref) which wraps:

1. A `PeriodicGraph`
2. The list of positions of the vertices in a unit cell of the graph
3. Optionally, a [`Cell`](@ref) if the dimension of the graph is 3 or below, which contains
   the geometry of the unit cell.

A `PeriodicGraphEmbedding` can be built through different methods, depending on whether the
list of positions should be permuted to be sorted, or offset to have all positions between
0 and 1 for instance:

```@docs
PeriodicGraphEmbedding
PeriodicGraphEmbedding3D
PeriodicGraphEmbedding{D,T}(graph::PeriodicGraph{D}, placement::AbstractMatrix{T}, cell::Cell) where {D,T}
sort_periodicgraphembedding!
PeriodicGraphEmbedding{D,T}(pge::PeriodicGraphEmbedding{N,S}) where {D,T,N,S}
```

## Cell API

```@docs
Cell
cell_parameters
EquivalentPosition
Base.parse(::Type{EquivalentPosition}, s::AbstractString)
find_refid
```
