# Main types of PeriodicGraphEmbeddings.jl

export EquivalentPosition,
       Cell,
       cell_parameters,
       SymmetryGroup3D,
       find_refid,
       PeriodicGraphEmbedding,
       PeriodicGraphEmbedding3D
       

## EquivalentPosition

"""
    EquivalentPosition

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
struct EquivalentPosition
    mat::SMatrix{3,3,Int,9}
    ofs::SVector{3,Rational{Int}}
end

(eq::EquivalentPosition)(x) = muladd(eq.mat, x, eq.ofs)

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

Parse a string into its represented `EquivalentPosition` given the name of the three
variables obtained from [`find_refid`](@ref PeriodicGraphEmbeddings.find_refid).

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
    something_written = false
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
                something_written = true
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
                    something_written || error(lazy"{$s} is not a symmetry equivalent (no dependency expressed in position $i)")
                    something_written = false
                    i += 1
                else
                    k !== Tokenize.Tokens.ENDMARKER && error(lazy"Unknown end of line marker for symmetry equivalent {$s}")
                    i != 3 && error(lazy"Input string \"$s\" is not a valid symmetry equivalent")
                end
            end
        end
    end
    EquivalentPosition(SMatrix{3,3,Int,9}(mat), SVector{3,Rational{Int}}(ofs))
end

function rationaltostring(x::Union{Int,Rational}, notofs::Bool, first::Bool)
    if notofs && (x == 1 || x == -1)
        return x < 0 ? "-" : first ? "" : "+"
    end
    nopluscond = x < 0 || first
    x isa Int && return nopluscond ? string(x) : string('+', x)
    num = numerator(x)
    den = denominator(x)
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
                coeff = rationaltostring(eq.mat[i,j], true, first)
                first = false
                print(io, coeff)
                print(io, xyz[j])
            end
        end
        if eq.ofs[i] != 0
            print(io, rationaltostring(eq.ofs[i], false, first))
        end
        i < 3 && print(io, ',')
    end
    nothing
end


## Cell

"""
    Cell

Representation of a periodic cell in 3D. Contains information about the cell
(axes lengths and angles) and its symmetry group, through its Hall number.

See [`SPACE_GROUP_HALL`](@ref),
[`SPACE_GROUP_FULL`](@ref),
[`SPACE_GROUP_HM`](@ref)
and [`SPACE_GROUP_IT`](@ref)
for the correspondance between Hall number and usual symbolic representations.
"""
struct Cell
    hall::Int
    mat::SMatrix{3,3,BigFloat,9} # Cartesian coordinates of a, b and c
    equivalents::Vector{EquivalentPosition}

    Cell(hall, mat::AbstractMatrix, equivalents) = new(hall, mat, equivalents)
end

function ==(c1::Cell, c2::Cell)
    c1.hall == c2.hall && c1.mat == c2.mat
end

function Cell(hall, (a, b, c), (α, β, γ), eqs=EquivalentPosition[])
    cosα = cosd(α); cosβ = cosd(β); cosγ = cosd(γ); sinγ = sind(γ)
    ω = sqrt(1 - cosα^2 - cosβ^2 - cosγ^2 + 2*cosα*cosβ*cosγ)
    mat = SMatrix{3,3,BigFloat,9}([a  b*cosγ  c*cosβ ;
                                  0   b*sinγ  c*(cosα - cosβ*cosγ)/sinγ ;
                                  0   0       c*ω/sinγ ])
    if isempty(eqs)
        eqs = get_symmetry_equivalents(hall)
        popfirst!(eqs) # corresponds to the identity
    end
    return Cell(hall, mat, eqs)
end
Cell() = Cell(1, (10, 10, 10), (90, 90, 90), EquivalentPosition[])

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

Cell(cell::Cell, mat::StaticArray{Tuple{3,3},BigFloat}) = Cell(cell.hall, mat, cell.equivalents)
Cell(mat::StaticArray{Tuple{3,3},BigFloat}) = Cell(Cell(), mat)

function Base.show(io::IO, cell::Cell)
    ((__a, __b, __c), (__α, __β, __γ)), _ = cell_parameters(cell)
    _a, _b, _c, _α, _β, _γ = Float64.((__a, __b, __c, __α, __β, __γ))
    print(io, Cell, '(', cell.hall, ", (", _a, ", ", _b, ", ", _c, "), (", _α, ", ", _β, ", ", _γ, ')')
end
function Base.show(io::IO, ::MIME"text/plain", cell::Cell)
    ((__a, __b, __c), (__α, __β, __γ)), _ = cell_parameters(cell)
    _a, _b, _c, _α, _β, _γ = Float64.((__a, __b, __c, __α, __β, __γ))
    hall_symbol, crystal_system = HALL_SYMBOLS[cell.hall]
    print(io, Cell, " with Hall symbol ", hall_symbol, " (", crystal_system, ") and parameters ")
    print(io, "a=", _a, ", b=", _b, ", c=", _c, ", α=", _α, ", β=", _β, ", γ=", _γ)
end


## PeriodicGraphEmbedding

"""
    PeriodicGraphEmbedding{D,T}

Embedding in euclidean space of a `PeriodicGraph` of dimension `D`.
Each vertex is assigned a `D`-uplet of coordinates of type `T`.

`PeriodicGraphEmbedding3D` is provided as an alias for `PeriodicGraphEmbedding{3}`.
Symmetry detection can only be performed on `PeriodicGraphEmbedding3D`.
"""
struct PeriodicGraphEmbedding{D,T}
    g::PeriodicGraph{D}
    pos::Vector{SVector{D,T}}
    cell::Cell

    function PeriodicGraphEmbedding{D,T}(g::PeriodicGraph{D}, pos::Vector{SVector{D,T}}, cell::Cell) where {D,T}
        new{D,T}(g, pos, cell)
    end
end

const PeriodicGraphEmbedding3D = PeriodicGraphEmbedding{3}

function PeriodicGraphEmbedding{D,T}(g::PeriodicGraph{D}, pos, cell::Cell=Cell()) where {D,T}
    if eltype(pos) == SVector{D,T}
        return PeriodicGraphEmbedding{D,T}(g, collect(pos), cell)
    end
    PeriodicGraphEmbedding{D,T}(g, [SVector{D,T}(p) for p in pos], cell)
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

function PeriodicGraphEmbedding{D,T}(graph::PeriodicGraph{D}, placement::AbstractMatrix{T}, cell::Cell, ::Val{order}=Val(false)) where {D,T,order}
    n = nv(graph)
    pos = Vector{SVector{D,T}}(undef, n)
    offsets = Vector{SVector{D,Int}}(undef, n)
    @inbounds for (i, x) in enumerate(eachcol(placement))
        offsets[i] = floor.(Int, x)
        pos[i] = Base.unsafe_rational.(getfield.(x, :num) .- getfield.(x, :den).*offsets[i], getfield.(x, :den))
    end
    if order
        s = sortperm(pos)
        pos = pos[s]
        graph = offset_representatives!(graph, .-offsets)[s]
        return PeriodicGraphEmbedding{D,T}(graph, pos, cell), s
    else
        graph = offset_representatives!(graph, .-offsets)
        return PeriodicGraphEmbedding{D,T}(graph, pos, cell)
    end
end

Base.getindex(pge::PeriodicGraphEmbedding, i::Integer) = pge.pos[i]
Base.getindex(pge::PeriodicGraphEmbedding, u::PeriodicVertex) = pge.pos[u.v] .+ u.ofs

function Base.getindex(pge::PeriodicGraphEmbedding{D,T}, vmap) where {D,T}
    PeriodicGraphEmbedding{D,T}(pge.g[vmap], pge.pos[vmap], pge.cell)
end


## SymmetryGroup3D

struct PeriodicSymmetry3D{T} <: PeriodicGraphs.AbstractSymmetry
    vmap::SubArray{PeriodicVertex3D,1,Matrix{PeriodicVertex3D},Tuple{Base.Slice{Base.OneTo{Int}},Int},true}
    rotation::SMatrix{3,3,Int,9}
    translation::SVector{3,T}
end
Base.getindex(symm::PeriodicSymmetry3D, i::Integer) = symm.vmap[i].v
function Base.getindex(symm::PeriodicSymmetry3D, x::PeriodicVertex3D)
    dst = symm.vmap[x.v]
    _ofs = muladd(symm.rotation, x.ofs, dst.ofs)
    PeriodicVertex3D(dst.v, _ofs)
end


"""
    SymmetryGroup3D{T} <: PeriodicGraphs.AbstractSymmetryGroup

Store the information on the available symmetry operations available on a
`PeriodicGraphEmbedding3D`.
"""
struct SymmetryGroup3D{T} <: PeriodicGraphs.AbstractSymmetryGroup{PeriodicSymmetry3D{T}}
    vmaps::Matrix{PeriodicVertex3D}
    rotations::Vector{SMatrix{3,3,Int,9}}
    translations::Vector{SVector{3,T}}
    uniquemap::Vector{Int}
    uniques::Vector{Int}
    hasmirror::Bool
end

(s::SymmetryGroup3D)(i::Integer) = s.uniques[s.uniquemap[i]]
Base.unique(s::SymmetryGroup3D) = s.uniques
function Base.getindex(s::SymmetryGroup3D{T}, i::Integer) where {T}
    PeriodicSymmetry3D{T}((@view s.vmaps[:,i]), s.rotations[i], s.translations[i])
end
Base.iterate(s::SymmetryGroup3D, state=1) = state > length(s) ? nothing : (s[state], state+1)
Base.eltype(::Type{SymmetryGroup3D}) = SubArray{Int,1,Matrix{Int},Tuple{Base.Slice{Base.OneTo{Int}},Int},true}
Base.length(s::SymmetryGroup3D) = size(s.vmaps, 2)
function Base.one(s::SymmetryGroup3D)
    n = length(s.uniquemap)
    @view reshape(collect(Base.OneTo(n)), n, 1)[:,1]
end

function SymmetryGroup3D{T}(vmaps_list, rots, transs, hasmirror, n) where T
    m = length(rots)
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
    return SymmetryGroup3D{T}(vmaps, rots, transs, uniquemap, uniques, hasmirror)
end
