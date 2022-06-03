# Main types of PeriodicGraphEmbeddings.jl

export EquivalentPosition,
       Cell,
       cell_parameters,
       SymmetryGroup3D,
       find_refid,
       PeriodicGraphEmbedding,
       PeriodicGraphEmbedding3D,
       SortedPeriodicGraphEmbedding


## EquivalentPosition

"""
    EquivalentPosition{T}

Representation of a symmetry operation in 3D, defined by a matrix multiplication and addition.

## Example
```jldoctest
julia> eq = parse(EquivalentPosition, "1-x, z, y+1/2")
-x+1,z,y+1/2

julia> eq([1//3, 0, 1//4])
3-element StaticArrays.SVector{3, Rational{$Int}} with indices SOneTo(3):
 2//3
 1//4
 1//2
```
"""
struct EquivalentPosition{T}
    mat::SMatrix{3,3,Int,9}
    ofs::SVector{3,T}
end

(eq::EquivalentPosition)(x) = muladd(eq.mat, x, eq.ofs)
EquivalentPosition{T}() where T = EquivalentPosition{T}(one(SMatrix{3,3,Int,9}), zero(SVector{3,T}))

"""
    find_refid(eqs)

Find the reference identifiers for the three dimensions for the CIF group called
`symmetry_equiv_pos_as_xyz` or `space_group_symop_operation_xyz`.
Usually this is simply ("x", "y", "z").
"""
function find_refid(eqs)
    isempty(eqs) && return ("x", "y", "z")
    for eq in eqs # Normally it should be the first one but I'm not sure of the specs here
        ts = collect(tokenize(eq))
        refid = [""]
        not_at_the_end = true
        immediate_continue = false
        for t in ts
            k = Tokenize.Tokens.kind(t)
            if k === Tokenize.Tokens.IDENTIFIER
                refid[end] = Tokenize.Tokens.untokenize(t)
            elseif k === Tokenize.Tokens.COMMA || k === Tokenize.Tokens.SEMICOLON
                push!(refid, "")
            elseif k === Tokenize.Tokens.ENDMARKER
                not_at_the_end = false
            elseif k !== Tokenize.Tokens.WHITESPACE && t.kind !== Tokenize.Tokens.PLUS
                immediate_continue = true
                break
            end
        end
        immediate_continue && continue
        if not_at_the_end
            error(lazy"Unknown end of line marker for symmetry equivalent {$eq}")
        end
        if length(refid) != 3 || refid[end] == ""
            error(lazy"Input string {$eq} is not a valid symmetry equivalent")
        end
        return tuple(lowercase(refid[1]), lowercase(refid[2]), lowercase(refid[3]))
    end
    return ("x", "y", "z")
end

"""
    Base.parse(::Type{EquivalentPosition}, s::AbstractString, refid=("x", "y", "z"))

Parse a string into its represented [`EquivalentPosition`](@ref) given the name of the
three variables obtained from [`find_refid`](@ref PeriodicGraphEmbeddings.find_refid).

## Example
```jldoctest
julia> parse(EquivalentPosition, "a+0.5; c; -1/3+b", ("a", "b", "c"))
x+1/2,z,y-1/3

julia> 
```
"""
function Base.parse(::Type{EquivalentPosition}, s::AbstractString, refid=("x", "y", "z"))
    const_dict = Dict{String, Int}(refid[1]=>1, refid[2]=>2, refid[3]=>3)
    mat = zeros(Int, 3, 3)
    ofs = zeros(Rational{Int}, 3)
    curr_num::Union{Int, Nothing} = nothing
    curr_val::Union{Rational{Int}, Nothing} = nothing
    curr_sign::Union{Bool, Nothing} = nothing
    encountered_div::Bool = false
    i::Int = 1
    for x in tokenize(lowercase(s))
        k = Tokenize.Tokens.kind(x)
        k === Tokenize.Tokens.WHITESPACE && continue
        if k === Tokenize.Tokens.INTEGER
            if encountered_div
                curr_val = Rational{Int}(Int(curr_num), parse(Int, x.val))
                curr_num = nothing
                encountered_div = false
            else
                curr_num = parse(Int, x.val)
            end
        else
            if k === Tokenize.Tokens.IDENTIFIER
                if !isnothing(curr_num)
                    curr_val = curr_num
                    curr_num = nothing
                end
                sign = isnothing(curr_sign) ? 1 : 2*curr_sign - 1
                val = isnothing(curr_val)  ? 1//1 : curr_val
                j = const_dict[Tokenize.Tokens.untokenize(x)]
                mat[i,j] += sign * val
                curr_val = nothing
                curr_sign = nothing
            else
                if x.kind === Tokenize.Tokens.FWD_SLASH
                    encountered_div = true
                    continue
                end
                if x.kind === Tokenize.Tokens.FLOAT
                    curr_val = Rational{Int}(rationalize(Int8, parse(Float16, x.val)))
                    continue
                end
                if !isnothing(curr_num)
                    curr_val = curr_num
                    curr_num = nothing
                end
                if !isnothing(curr_val)
                    sign = isnothing(curr_sign) ? 1 : 2*curr_sign - 1
                    ofs[i] += sign * Rational{Int}(curr_val)
                    curr_val = nothing
                    curr_sign = nothing
                end
                if x.kind === Tokenize.Tokens.PLUS
                    curr_sign = true
                elseif x.kind === Tokenize.Tokens.MINUS
                    curr_sign = false
                elseif k === Tokenize.Tokens.COMMA || k === Tokenize.Tokens.SEMICOLON
                    i > 2 && error(lazy"Too many dimensions specified for symmetry equivalent {$s}")
                    i += 1
                else
                    k !== Tokenize.Tokens.ENDMARKER && error(lazy"Unknown end of line marker for symmetry equivalent {$s}")
                    i != 3 && error(lazy"Input string \"$s\" is not a valid symmetry equivalent")
                end
            end
        end
    end
    EquivalentPosition{Rational{Int}}(SMatrix{3,3,Int,9}(mat), SVector{3,Rational{Int}}(ofs))
end

function __tostring(x::T, notofs::Bool, first::Bool) where T
    if notofs && (x == 1 || x == -1)
        return x < 0 ? "-" : first ? "" : "+"
    end
    nopluscond = x < 0 || first
    T == Int && return nopluscond ? string(x) : string('+', x)
    num, den = T <: Rational ? (numerator(x), denominator(x)) : (x, 1)
    if isone(den)
        return nopluscond ? string(num) : string('+', num)
    end
    return nopluscond ? string(num, '/', den) : string('+', num, '/', den)
end

function Base.show(io::IO, eq::EquivalentPosition)
    xyz = ('x', 'y', 'z')
    for i in 1:3
        first = true
        for j in 1:3
            if eq.mat[i,j] != 0
                coeff = __tostring(eq.mat[i,j], true, first)
                first = false
                print(io, coeff)
                print(io, xyz[j])
            end
        end
        if eq.ofs[i] != 0 || first
            print(io, __tostring(eq.ofs[i], false, first))
        end
        i < 3 && print(io, ',')
    end
    nothing
end


## Cell

"""
    Cell{T}

Representation of a periodic cell in 3D. Contains information about the cell
(axes lengths and angles) and its symmetry group, through its Hall number.

See [`PeriodicGraphEmbeddings.SPACE_GROUP_HALL`](@ref),
[`PeriodicGraphEmbeddings.SPACE_GROUP_FULL`](@ref),
[`PeriodicGraphEmbeddings.SPACE_GROUP_HM`](@ref)
and [`PeriodicGraphEmbeddings.SPACE_GROUP_IT`](@ref)
for the correspondance between Hall number and usual symbolic representations.
"""
struct Cell{T}
    hall::Int
    mat::SMatrix{3,3,BigFloat,9} # Cartesian coordinates of a, b and c
    equivalents::Vector{EquivalentPosition{T}}

    Cell{T}(hall, mat::AbstractMatrix, equivalents) where {T} = new{T}(hall, mat, equivalents)
    function Cell(hall, mat::AbstractMatrix, equivalents::AbstractVector{EquivalentPosition{T}}) where T
        new{T}(hall, mat, equivalents)
    end
end

function Cell{T}(hall, (a, b, c), (α, β, γ), eqs=EquivalentPosition{T}[]) where T
    cosα = cosd(α); cosβ = cosd(β); cosγ = cosd(γ); sinγ = sind(γ)
    ω = sqrt(1 - cosα^2 - cosβ^2 - cosγ^2 + 2*cosα*cosβ*cosγ)
    mat = SMatrix{3,3,BigFloat,9}([a  b*cosγ  c*cosβ ;
                                  0   b*sinγ  c*(cosα - cosβ*cosγ)/sinγ ;
                                  0   0       c*ω/sinγ ])
    if isempty(eqs)
        eqs = get_symmetry_equivalents(T, hall)
        popfirst!(eqs) # corresponds to the identity
    end
    return Cell{T}(hall, mat, eqs)
end
Cell(hall, vecs, angs, eqs=EquivalentPosition{Rational{Int}}[]) = Cell{Rational{Int}}(hall, vecs, angs, eqs)
Cell{T}() where {T} = Cell(1, (10, 10, 10), (90, 90, 90), EquivalentPosition{T}[])
Cell() = Cell{Rational{Int}}()

Cell(cell::Cell{T}, mat::StaticArray{Tuple{3,3},BigFloat}) where {T} = Cell{T}(cell.hall, mat, cell.equivalents)
Cell{T}(mat::StaticArray{Tuple{3,3},BigFloat}) where {T} = Cell(Cell{T}(), mat)
Cell(mat::StaticArray{Tuple{3,3},BigFloat}) = Cell{Rational{Int}}(mat)

function cell_parameters(mat::AbstractMatrix)
    _a, _b, _c = eachcol(mat)
    a = norm(_a)
    b = norm(_b)
    c = norm(_c)
    α = acosd(_b'_c/(b*c))
    β = acosd(_c'_a/(c*a))
    γ = acosd(_a'_b/(a*b))
    return (a, b, c), (α, β, γ)
end

"""
    cell_parameters(cell::Cell)

Return `((lengths, angles), mat)` where `mat` is the matrix of the cell in upper triangular
format, `lengths` is the triplet `(a, b, c)` of lengths of the three axes, and `angles` is
the triplet `(α, β, γ)` of angles between them.
"""
function cell_parameters(cell::Cell)
    lengths, angles = cell_parameters(cell.mat)
    normalized_mat = istriu(cell.mat) ? cell.mat :
        Cell(cell.hall, lengths, angles, cell.equivalents).mat
    return (lengths, angles), normalized_mat
end

function ==(cell1::Cell, cell2::Cell)
    cell1.hall == cell2.hall || return false
    cell1.mat == cell2.mat && return true
    istriu(cell1.mat) && istriu(cell2.mat) && return false
    (a1, b1, c1), (α1, β1, γ1) = cell_parameters(cell1.mat)
    (a2, b2, c2), (α2, β2, γ2) = cell_parameters(cell2.mat)
    return Float32.((a1, b1, c1, α1, β1, γ1)) == Float32.((a2, b2, c2, α2, β2, γ2))
end


function Base.show(io::IO, cell::Cell)
    ((__a, __b, __c), (__α, __β, __γ)), _ = cell_parameters(cell)
    _a, _b, _c, _α, _β, _γ = Float64.((__a, __b, __c, __α, __β, __γ))
    print(io, cell isa Cell{Rational{Int}} ? Cell : typeof(cell))
    print(io, '(', cell.hall, ", (", _a, ", ", _b, ", ", _c, "), (", _α, ", ", _β, ", ", _γ, ')', ')')
end
function Base.show(io::IO, ::MIME"text/plain", cell::Cell)
    ((__a, __b, __c), (__α, __β, __γ)), _ = cell_parameters(cell)
    _a, _b, _c, _α, _β, _γ = Float64.((__a, __b, __c, __α, __β, __γ))
    hall_symbol, crystal_system = HALL_SYMBOLS[cell.hall]
    print(io, cell isa Cell{Rational{Int}} ? Cell : typeof(cell))
    print(io, " with Hall symbol ", hall_symbol, " (", crystal_system, ") and parameters ")
    print(io, "a=", _a, ", b=", _b, ", c=", _c, ", α=", _α, ", β=", _β, ", γ=", _γ)
end


## PeriodicGraphEmbedding

"""
    PeriodicGraphEmbedding{D,T}

Embedding in euclidean space of a `PeriodicGraph` of dimension `D`.
Each vertex is assigned a `D`-uplet of coordinates of type `T`.

`PeriodicGraphEmbedding3D` is provided as an alias for `PeriodicGraphEmbedding{3}`.
Symmetry detection provided by PeriodicGraphEmbeddings.jl can only be performed on
`PeriodicGraphEmbedding3D`.
"""
struct PeriodicGraphEmbedding{D,T}
    g::PeriodicGraph{D}
    pos::Vector{SVector{D,T}}
    cell::Cell{Rational{Int}}

    function PeriodicGraphEmbedding{D,T}(g::PeriodicGraph{D}, pos::Vector{SVector{D,T}}, cell::Cell{Rational{Int}}) where {D,T}
        new{D,T}(g, pos, cell)
    end
end

==(x::PeriodicGraphEmbedding, y::PeriodicGraphEmbedding) = x.g == y.g && x.pos == y.pos && x.cell == y.cell

"""
    PeriodicGraphEmbedding3D

Alias for `PeriodicGraphEmbedding{3}`
"""
const PeriodicGraphEmbedding3D = PeriodicGraphEmbedding{3}

function PeriodicGraphEmbedding{D,T}(g::PeriodicGraph{D}, pos, cell::Cell=Cell()) where {D,T}
    if eltype(pos) == SVector{D,T}
        return PeriodicGraphEmbedding{D,T}(g, collect(pos)::Vector{SVector{D,T}}, cell)
    end
    PeriodicGraphEmbedding{D,T}(g, SVector{D,T}[SVector{D,T}(p) for p in pos], cell)
end
function PeriodicGraphEmbedding{D}(g, pos::S, cell::Cell=Cell()) where {D,S}
    elt = eltype(S)
    if Base.IteratorEltype(elt) isa Base.HasEltype
        return PeriodicGraphEmbedding{D,eltype(elt)}(g, pos, cell)
    else
        throw(MethodError(PeriodicGraphEmbedding{D}, typeof((g, pos, cell))))
    end
end
function PeriodicGraphEmbedding(g::PeriodicGraph{D}, pos, cell::Cell=Cell()) where D
    PeriodicGraphEmbedding{D}(g, pos, cell)
end

function PeriodicGraphEmbedding{D,T}(cell::Cell=Cell()) where {D,T}
    PeriodicGraphEmbedding{D,T}(PeriodicGraph{D}(), SVector{D,T}[], cell)
end

function return_to_box(g::PeriodicGraph{D}, placement::AbstractMatrix{T}) where {D,T}
    n = nv(g)
    pos = Vector{SVector{D,T}}(undef, n)
    offsets = Vector{SVector{D,Int}}(undef, n)
    @inbounds for (i, x) in enumerate(eachcol(placement))
        ofs = floor.(Int, x)
        offsets[i] = ofs
        pos[i] = @. Base.unsafe_rational(numerator(x) - denominator(x)*ofs, denominator(x))
    end
    return pos, offsets
end

"""
    PeriodicGraphEmbedding{D,T}(graph::PeriodicGraph{D}, placement::AbstractMatrix{T}, cell::Cell=Cell()) where {D,T}
    PeriodicGraphEmbedding{D}(graph::PeriodicGraph{D}, placement::AbstractMatrix{T}, cell::Cell=Cell()) where D
    PeriodicGraphEmbedding(graph::PeriodicGraph{D}, placement::AbstractMatrix{T}, cell::Cell=Cell())

Build a [`PeriodicGraphEmbedding{D,T}`](@ref) from the corresponding `graph` and
`placement` of the vertices, such that each vertex has its fractional coordinate
represented in a column of the matrix.

Coordinates out of [0, 1) are translated back to the unit cell with the corresponding
offset added to the graph.

The `cell` optional argument will not be used if `D > 3`.

!!! warning
    This function modifies the input `graph` if any element of `placement` is out of [0, 1).

!!! note
    To obtain a [`PeriodicGraphEmbedding`](@ref) with sorted positions, use
    [`SortedPeriodicGraphEmbedding`](@ref) instead
"""
function PeriodicGraphEmbedding{D,T}(graph::PeriodicGraph{D}, placement::AbstractMatrix{T}, cell::Cell) where {D,T}
    pos, offsets = return_to_box(graph, placement)
    offset_representatives!(graph, .-offsets)
    return PeriodicGraphEmbedding{D,T}(graph, pos, cell)
end

"""
    SortedPeriodicGraphEmbedding{T}

Constructor for `PeriodicGraphEmbedding{D,T} where D` with sorted positions.
"""
struct SortedPeriodicGraphEmbedding{T}
    __noconstructor() = nothing
end

"""
    SortedPeriodicGraphEmbedding{T}(graph::PeriodicGraph{D}, placement::AbstractMatrix, cell::Cell=Cell()) where {D,T}

Build a [`PeriodicGraphEmbedding{D,T}`](@ref) from the corresponding `graph` and
`placement` of the vertices, so that the result has its vertices sorted by position.

Return the [`PeriodicGraphEmbedding`](@ref) as well as the permutation of the columns of
`placement` that yielded the resulting order on the vertices.

The `cell` optional argument will not be used if `D > 3`.

!!! warning
    This function modifies the input `graph` if any element of `placement` is out of [0, 1).

See also [`PeriodicGraphEmbedding{D,T}(graph, placement::AbstractMatrix{T}, cell) where {D,T}`](@ref)
and [`SortedPeriodicGraphEmbedding(graph, placement::AbstractMatrix, cell)`](@ref).
"""
function SortedPeriodicGraphEmbedding{T}(graph::PeriodicGraph{D}, placement::AbstractMatrix, cell::Cell=Cell()) where {D,T}
    pos, offsets = return_to_box(graph, placement)
    s = sortperm(pos)
    if s == 1:(length(pos)) # avoid rebuilding the graph if unnecessary
        offset_representatives!(graph, .-offsets)
        return PeriodicGraphEmbedding{D,T}(graph, pos, cell), s
    end
    pos = pos[s]
    graph = offset_representatives!(graph, .-offsets)[s]
    return PeriodicGraphEmbedding{D,T}(graph, pos, cell), s
end

"""
    double_widen(::Type)

Internal function used to selectively widen small integer and rational types.

This is useful to avoid overflow without sacrificing too much efficiency by
always having to resolve to very large types.
"""
double_widen(::Type{T}) where {T} = T
double_widen(::Type{Int64}) = Int128
double_widen(::Type{Int32}) = Int64
double_widen(::Type{Int16}) = Int64
double_widen(::Type{Int8}) = Int32

macro tryinttype(T)
    tmin = :(typemin($T))
    tmax = :(typemax($T))
    S = :(double_widen($T))
    return esc(quote
        if (($tmin <= m) & (M <= $tmax))
            return SortedPeriodicGraphEmbedding{Rational{$S}}(graph, placement, cell)
        end
    end)
end

"""
    SortedPeriodicGraphEmbedding(graph::PeriodicGraph{D}, placement::AbstractMatrix, cell::Cell=Cell()) where D

Build a [`PeriodicGraphEmbedding{D,T}`](@ref) from the corresponding `graph` and `placement` of the
vertices, so that the result has its vertices sorted by position. `T` is determined as the
smallest type between `Rational{Int32}`, `Rational{Int64}`, `Rational{Int128}` and
`Rational{BigInt}` that can fit all the elements of `placement` with some additional
margin.

Return the [`PeriodicGraphEmbedding`](@ref) as well as the permutation of the columns of
`placement` that yielded the resulting order on the vertices.

The `cell` optional argument will not be used if `D > 3`.

!!! warning
    This function modifies the input `graph` if any element of `placement` is out of [0, 1).

!!! tip
    This function is inherently type-unstable since `T` cannot be statically determined.
    This can be useful because having a too large `T` may slow down later computations.
    
    To provide the parameter explicitly, pass it to the `SortedPeriodicGraphEmbedding`
    constructor by calling `SortedPeriodicGraphEmbedding{T}(graph, placement, cell)`.

See also [`PeriodicGraphEmbedding{D,T}(graph, placement::AbstractMatrix{T}, cell) where {D,T}`](@ref).
"""
function SortedPeriodicGraphEmbedding(graph::PeriodicGraph{D}, placement::AbstractMatrix, cell::Cell=Cell()) where D
    if isempty(placement)
        return PeriodicGraphEmbedding{D,Rational{Int32}}(cell), Int[]
    end
    m = min(minimum(numerator.(placement)), minimum(denominator.(placement)))
    M = max(maximum(numerator.(placement)), maximum(denominator.(placement)))
    @tryinttype Int8
    @tryinttype Int16
    @tryinttype Int32
    @tryinttype Int64
    @tryinttype Int128
    return SortedPeriodicGraphEmbedding{Rational{BigInt}}(graph, placement, cell)
    # Type-unstable function, but yields better performance than always falling back to BigInt
end


Base.getindex(pge::PeriodicGraphEmbedding, i::Integer) = pge.pos[i]
Base.getindex(pge::PeriodicGraphEmbedding, u::PeriodicVertex) = pge.pos[u.v] .+ u.ofs
Base.length(pge::PeriodicGraphEmbedding) = length(pge.pos)

function Base.getindex(pge::PeriodicGraphEmbedding{D,T}, vmap) where {D,T}
    PeriodicGraphEmbedding{D,T}(pge.g[vmap], pge.pos[vmap], pge.cell)
end

"""
    PeriodicGraphEmbedding{D,T}(pge::PeriodicGraphEmbedding{N,S}) where {D,T,N,S}
    PeriodicGraphEmbedding{D}(pge::PeriodicGraphEmbedding{N,S}) where {D,N,S}

Return a [`PeriodicGraphEmbedding{D,T}`](@ref) with the same structural information as the
input `pge` but embedded in `D` dimensions instead of `N`.

If `T` is not provided it defaults to `S`.

The same caveats that apply to `PeriodicGraph{D}(graph::PeriodicGraph{N})`
are valid here: namely, the dimensionality of the graph should be at least `D` and the
behaviour is undefined if `D < N` and there are multiple non-identical connected components.

Moreover, if `D < N`, the `N-D` last coordinates of all vertices must be
zero or this function will error.
"""
function PeriodicGraphEmbedding{D,T}(pge::PeriodicGraphEmbedding{N,S}) where {D,T,N,S}
    N == D && S == T && return pge
    newpositions = Vector{SVector{D,T}}(undef, length(pge))
    if D > N
        @simd for i in 1:length(pge)
            @inbounds newpositions[i] = SVector{D,T}([pge.pos[i]; zero(SVector{D-N,T})])
        end
    else
        for (i, pos) in enumerate(pge.pos)
            for j in (D+1):N
                pos[j] == 0 || throw(DimensionMismatch(lazy"Cannot convert to dimension D=$D because coordinate $j of vertex $i is not 0 and $j > D."))
            end
            @inbounds newpositions[i] = SVector{D,T}(@view pos[1:D])
        end
    end
    newg = PeriodicGraph{D}(pge.g)
    return PeriodicGraphEmbedding{D,T}(newg, newpositions, pge.cell)
end
PeriodicGraphEmbedding{D}(pge::PeriodicGraphEmbedding{N,S}) where {D,N,S} = PeriodicGraphEmbedding{D,S}(pge)


## SymmetryGroup3D

"""
    PeriodicSymmetry3D{T} <: PeriodicGraphs.AbstractSymmetry

Single symmetry of a [`PeriodicGraphEmbedding3D{T}`](@ref).

See [`PeriodicGraphs.AbstractSymmetry`](https://liozou.github.io/PeriodicGraphs.jl/dev/symmetries/#PeriodicGraphs.AbstractSymmetry)
for information on the API.
"""
struct PeriodicSymmetry3D{T} <: PeriodicGraphs.AbstractSymmetry
    vmap::SubArray{PeriodicVertex3D,1,Matrix{PeriodicVertex3D},Tuple{Base.Slice{Base.OneTo{Int}},Int},true}
    eq::EquivalentPosition{T}
end
Base.getindex(symm::PeriodicSymmetry3D, i::Integer) = symm.vmap[i].v
function Base.getindex(symm::PeriodicSymmetry3D, x::PeriodicVertex3D)
    dst = symm.vmap[x.v]
    _ofs = muladd(symm.eq.mat, x.ofs, dst.ofs)
    PeriodicVertex3D(dst.v, _ofs)
end
(symm::PeriodicSymmetry3D)(x::AbstractVector) = (symm.eq)(x)
(symm::PeriodicSymmetry3D)(x::AbstractMatrix) = symm.eq.mat * x


"""
    SymmetryGroup3D{T} <: PeriodicGraphs.AbstractSymmetryGroup

Store the information on the symmetry operations available on a
[`PeriodicGraphEmbedding3D`](@ref).
"""
struct SymmetryGroup3D{T} <: PeriodicGraphs.AbstractSymmetryGroup{PeriodicSymmetry3D{T}}
    vmaps::Matrix{PeriodicVertex3D}
    eqs::Vector{EquivalentPosition{T}}
    uniquemap::Vector{Int}
    uniques::Vector{Int}
    hasmirror::Bool
end

(s::SymmetryGroup3D)(i::Integer) = s.uniques[s.uniquemap[i]]
Base.unique(s::SymmetryGroup3D) = s.uniques
function Base.getindex(s::SymmetryGroup3D{T}, i::Integer) where {T}
    PeriodicSymmetry3D{T}((@view s.vmaps[:,i]), s.eqs[i])
end
Base.iterate(s::SymmetryGroup3D, state=1) = state > length(s) ? nothing : (s[state], state+1)
Base.length(s::SymmetryGroup3D) = size(s.vmaps, 2)
function Base.one(s::SymmetryGroup3D{T}) where T
    n = length(s.uniquemap)
    mat = reshape([PeriodicVertex3D(i) for i in 1:n], n, 1)
    PeriodicSymmetry3D{T}((@view mat[:,1]), EquivalentPosition{T}())
end

function SymmetryGroup3D(vmaps_list, eqs::Vector{EquivalentPosition{T}}, hasmirror, n) where T
    m = length(eqs)
    keepuniques = trues(n)
    vmaps = Matrix{PeriodicVertex3D}(undef, n, m)
    for i in 1:m
        listi = vmaps_list[i]
        vmaps[:,i] .= vmaps_list[i]
        for j in 2:n
            keepuniques[j] &= listi[j].v ≥ j
        end
    end
    uniques = (1:n)[keepuniques]
    uniquemap = Vector{Int}(undef, n)
    for (k,j) in enumerate(uniques)
        uniquemap[j] = k
    end
    undetermined = [i for i in 2:n if !keepuniques[i]]
    undeterminedmask = Vector{Int}(undef, length(undetermined)÷2) # sizehint
    for vmap in eachcol(vmaps)
        empty!(undeterminedmask)
        for (k,j) in enumerate(undetermined)
            vmapj = vmap[j].v
            if keepuniques[vmapj]
                uniquemap[j] = uniquemap[vmapj]
            else
                push!(undeterminedmask, k)
            end
        end
        isempty(undeterminedmask) && break
        undetermined = undetermined[undeterminedmask]
    end
    return SymmetryGroup3D{T}(vmaps, eqs, uniquemap, uniques, hasmirror)
end

function SymmetryGroup3D{T}(vmaps_list, rots::AbstractVector{<:SMatrix}, transs, hasmirror, n) where T
    eqs = [EquivalentPosition{T}(r, t) for (r, t) in zip(rots, transs)]
    return SymmetryGroup3D(vmaps_list, eqs, hasmirror, n)
end
