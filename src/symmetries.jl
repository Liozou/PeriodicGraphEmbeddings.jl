# Wrapper around the spglib library

include("symmetry_groups.jl")

export find_hall_number,
       get_symmetry_equivalents,
       check_valid_symmetry,
       find_symmetries

import spglib_jll: libsymspg
import LinearAlgebra

"""
    SpglibDataset

Wrapper around the `SpglibDataset` type exported by spglib.
Its accessible fields are the same as in the C counterpart, except that strings are already
converted to `String`, lists to `Vector` and matrices to `Matrix`.

To access the raw pointers without conversion, prepend an underscore to the field: for
example `dataset._rotations` yields a `Ptr{Cint}` where `dataset.rotations` is a 3×3
`Matrix{Int}`.
"""
mutable struct SpglibDataset
    spacegroup_number::Cint
    hall_number::Cint
    _international_symbol::NTuple{11,Cchar}
    _hall_symbol::NTuple{17,Cchar}
    _choice::NTuple{6,Cchar}
    _transformation_matrix::NTuple{9,Cdouble}
    _origin_shift::NTuple{3,Cdouble}
    n_operations::Cint
    _rotations::Ptr{Cint}
    _translations::Ptr{Cdouble}
    n_atoms::Cint
    _wyckoffs::Ptr{Cint}
    _site_symmetry_symbols::Ptr{NTuple{7,Cchar}}
    _equivalent_atoms::Ptr{Cint}
    _crystallographic_orbits::Ptr{Cint}
    _primitive_lattice::NTuple{9,Cdouble}
    _mapping_to_primitive::Ptr{Cint}
    n_std_atoms::Cint
    _std_lattice::NTuple{9,Cdouble}
    _std_types::Ptr{Cint}
    _std_positions::Ptr{Cdouble}
    _std_rotation_matrix::NTuple{9,Cdouble}
    _std_mapping_to_primitive::Ptr{Cint}
    _pointgroup_symbol::NTuple{6,Cchar}
end

function joinletters(x)
    letters = Char.(x)
    fzero = findfirst(==('\0'), letters)
    return join(letters[i] for i in 1:(isnothing(fzero) ? length(letters) : (fzero-1)))
end

function Base.getproperty(ds::SpglibDataset, name::Symbol)
    if name === :international_symbol
        return joinletters(getfield(ds, :_international_symbol))
    elseif name === :hall_symbol
        return joinletters(getfield(ds, :_hall_symbol))
    elseif name === :choice
        return joinletters(getfield(ds, :_choice))
    elseif name === :transformation_matrix
        return reshape(collect(getfield(ds, :_transformation_matrix)), (3,3))
    elseif name === :origin_shift
        return reshape(collect(getfield(ds, :_origin_shift)), 3)
    elseif name === :rotations
        return unsafe_wrap(Array, getfield(ds, :_rotations), (3,3,Int(getfield(ds, :n_operations))))
    elseif name === :translations
        return unsafe_wrap(Array, getfield(ds, :_translations), (3,Int(getfield(ds, :n_operations))))
    elseif name === :wyckoffs
        return unsafe_wrap(Array, getfield(ds, :_wyckoffs), Int(getfield(ds, :n_atoms)))
    elseif name === :site_symmetry_symbols
        a = unsafe_wrap(Array, getfield(ds, :_site_symmetry_symbols), getfield(ds, :n_atoms))
        ret = Vector{String}(undef, length(a))
        for i in 1:length(a)
            ret[i] = joinletters(a[i])
        end
        return ret
    elseif name === :equivalent_atoms
        return unsafe_wrap(Array, getfield(ds, :_equivalent_atoms), getfield(ds, :n_atoms))
    elseif name === :crystallographic_orbits
        return unsafe_wrap(Array, getfield(ds, :_crystallographic_orbits), getfield(ds, :n_atoms))
    elseif name === :primitive_lattice
        return reshape(collect(getfield(ds, :_primitive_lattice)), (3,3))
    elseif name === :mapping_to_primitive
        return unsafe_wrap(Array, getfield(ds, :_mapping_to_primitive), getfield(ds, :n_atoms))
    elseif name === :std_lattice
        return reshape(collect(getfield(ds, :_std_lattice)), (3,3))
    elseif name === :std_types
        return unsafe_wrap(Array, getfield(ds, :_std_types), getfield(ds, :n_std_atoms))
    elseif name === :std_positions
        return unsafe_wrap(Array, getfield(ds, :_std_positions), (3,Int(getfield(ds, :n_std_atoms))))
    elseif name === :std_rotation_matrix
        return reshape(collect(getfield(ds, :_std_rotation_matrix)), (3,3))
    elseif name === :std_mapping_to_primitive
        return unsafe_wrap(Array, getfield(ds, :_std_mapping_to_primitive), getfield(ds, :n_std_atoms))
    elseif name === :pointgroup_symbol
        return joinletters(getfield(ds, :_pointgroup_symbol))
    else
        getfield(ds, name)
    end
end


function find_hall_number(hallsymbol, hm=hallsymbol, it=0, warnonnotfound=false)
    if hallsymbol != ""
        hsymbol = lowercase(hallsymbol)
        hall = get(SPACE_GROUP_HALL, replace(hsymbol, ('_' => ' ')), 0)
        if hall ∈ (322, 326, 330)
            if hall == 322
                if endswith(hm, 'b')
                    hall = 324
                end
            elseif hall == 326
                if occursin('c', hm)
                    hall = 328
                end
            elseif occursin('a', hm)
                hall = 332
            end
        end
        hall == 0 || return hall
        warnonnotfound && @warn lazy"Hall symbol provided but not recognised: $hallsymbol"
    end
    if hm != ""
        hm = hm[1]*lowercase(@view hm[2:end])
        dense_hm = replace(join(split(hm)), ('_'=>""))
        hall = get(SPACE_GROUP_HM, dense_hm, 0)
        hall == 0 || return hall
        hall = get(SPACE_GROUP_FULL, dense_hm, 0)
        hall == 0 || return hall
        warnonnotfound && begin
            @warn lazy"H-M symbol provided but not recognised : $hm"
            keysHM = collect(keys(SPACE_GROUP_HM))
            suggestions = keysHM[findall(startswith(dense_hm), keysHM)]
            keysFULL = collect(keys(SPACE_GROUP_FULL))
            append!(suggestions, keysFULL[findall(startswith(dense_hm), keysFULL)])
            if !isempty(suggestions)
                msg = join(suggestions, ", ", " or ")
                @info lazy"Did you perhaps mean to use $msg?"
            end
        end
    end
    if it != 0
        if it < 1 || it > 230
            warnonnotfound && @error lazy"International Table number provided is outside of the allowed range (1-230) : $it"
            return 1
        end
        return SPACE_GROUP_IT[it]
    end
    warnonnotfound && @error "Could not determine the space group of this file. It will be considered P1."
    return 1
end



"""
    get_symmetry_equivalents(hall)

The list of `EquivalentPosition` corresponding to a symmetry group given by its Hall number.

Wrapper around `spg_get_symmetry_from_database`.
"""
function get_symmetry_equivalents(hall)
    rotations = Array{Cint}(undef, 3, 3, 192)
    translations = Array{Cdouble}(undef, 3, 192)
    len = ccall((:spg_get_symmetry_from_database, libsymspg), Cint,
                (Ptr{Cint}, Ptr{Cdouble}, Cint), rotations, translations, hall)
    eqs = EquivalentPosition[]
    for i in 1:len
        rot = SMatrix{3,3,Int,9}(transpose(@view rotations[:,:,i]))
        tr = SVector{3,Cdouble}(@view translations[:,i])
        trans = SVector{3,Rational{Int}}(round.(Int, 360 .* tr) .// 360)
        push!(eqs, EquivalentPosition(rot, trans))
    end
    return eqs
end

"""
    get_spglib_dataset(pge::PeriodicGraphEmbedding{3}, vtypes=nothing)

Wrapper around `spg_get_dataset`.

If `vtypes !== nothing`, ensure that two vertices `x` and `y` cannot be symmetry-related
if `vtypes[x] != vtypes[y]`.
"""
function get_spglib_dataset(pge::PeriodicGraphEmbedding{3,T}, vtypes=nothing) where T
    lattice = Matrix{Cdouble}(adjoint(pge.cell.mat)) # transpose to account for row-major operations
    n = nv(pge.graph)
    positions = Matrix{Cdouble}(undef, 3, n)
    if vtypes !== nothing
        types = Vector{Cint}(undef, n)
        vtypes_to_int = Dict{valtype(vtypes),Cint}()
        j = 1
    else
        types = zeros(Cint, n)
    end
    for i in 1:n
        positions[:,i] .= pge.cell.mat * pge.pos[i]
        if vtypes !== nothing
            j += ((types[i] = get!(vtypes_to_int, vtypes[i], j)) == j)
        end
    end
    ptr = ccall((:spg_get_dataset, libsymspg), Ptr{SpglibDataset},
                (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
                lattice, positions, types, n, T <: Rational ? 10*eps(Cdouble) : 0.0001)
    dataset = unsafe_load(ptr)
    # if dataset.n_atoms != n
    #     error("The periodic graph is not minimal") # TODO: update pge?
    # end
    return dataset
end


"""
    check_valid_symmetry(pge::PeriodicGraphEmbedding{D,T}, t::SVector{D,T}, r=nothing, vtypes=nothing, issorted=false)

Check that the periodic graph embedding is identical to that rotated by `r` (if it is not
`nothing`) then translated by `t`.
If `vtypes` is not `nothing`, any vertex `x` must additionally be mapped to a vertex `y`
such that `vtypes[x] == vtypes[y].`
If `issorted` is set and `T <: Rational`, assume that `issorted(pge.pos)` to use a faster
dichotomy approach.

If so, return the the `vmap` between the initial vertices and their symmetric images, as
well as the `offsets` of each symmetric image compared to the origin.
Otherwise, return `nothing`.
"""
function check_valid_symmetry(pge::PeriodicGraphEmbedding{D,T}, t::SVector{D,T}, r=nothing, vtypes=nothing, issorted=false) where {D,T}
    n = length(pge.pos)
    vmap = Vector{PeriodicVertex{D}}(undef, n)
    dichotomy = T <: Rational && issorted
    if !dichotomy
        notencountered = trues(n)
        if !(T <: Rational)
            mat = Float64.(pge.cell.mat)
            buffer, ortho, safemin = prepare_periodic_distance_computations(mat)
        end
    end
    for k in 1:n
        curr_pos = pge.pos[k]
        transl = (isnothing(r) ? curr_pos : (r * curr_pos)) .+ t
        ofs = floor.(Int, transl)
        x = transl .- ofs
        if dichotomy
            i, j = x < curr_pos ? (1, k) : (k, length(pge.pos)+1)
            while j - i > 1
                m = (j+i)>>1
                if cmp(x, pge.pos[m]) < 0
                    j = m
                else
                    i = m
                end
            end
            i -= i > length(pge.pos)
            pge.pos[i] == x || return nothing
        else
            i::Int = findfirst(notencountered)
            _j = findnext(notencountered, i+1)
            mindist = T <: Rational ? norm(pge.pos[i] .- x) : periodic_distance!(buffer, pge.pos[i] .- x, mat, ortho, safemin)
            while _j isa Int
                newdiff = pge.pos[_j] .- x
                newdist = T <: Rational ? norm(newdiff) : periodic_distance!(buffer, newdiff, mat, ortho, safemin)
                if newdist < mindist
                    mindist = newdist
                    i = _j
                end
                _j = findnext(notencountered, _j+1)
            end
            (T <: Rational ? pge.pos[i] == x : isapprox(mindist, 0.0; atol=0.0001)) || return nothing
            notencountered[i] = false
        end
        vtypes === nothing || vtypes[i] == vtypes[k] || return nothing
        vmap[k] = PeriodicVertex{D}(i, ofs)
    end
    for e in edges(pge.g)
        src = vmap[e.src]
        dst = vmap[e.dst.v]
        newofs = (isnothing(r) ? e.dst.ofs : r * e.dst.ofs) .+ dst.ofs .- src.ofs
        has_edge(pge.g, PeriodicGraphs.unsafe_edge{D}(src.v, dst.v, newofs)) || return nothing
    end
    return vmap
end

__rattype(::Type{Rational{T}}) where {T}  = T

"""
    find_symmetries(pge::PeriodicGraphEmbedding3D, vtypes=nothing, check_symmetry=check_valid_symmetry)

Return a `SymmetryGroup3D` object storing the list of symmetry operations on the graph embedding.

If `vtypes !== nothing`, ensure that two vertices `x` and `y` cannot be symmetry-related
if `vtypes[x] != vtypes[y]`.

`check_symmetry` must be a function that takes the same four arguments `pge`, `t`, `r` and
`vtypes` as `check_valid_symmetry` and return either `(vmap, offsets)` or `nothing` if the
input is not a valid symmetry.
It can be used to specify additional constraints that cannot be carried by `vtypes` alone.
"""
function find_symmetries(pge::PeriodicGraphEmbedding3D{T}, vtypes=nothing, check_symmetry=check_valid_symmetry) where T
    lattice = Matrix{Cdouble}(LinearAlgebra.I, 3, 3) # positions are expressed in this basis

    I = sortperm(pge.pos)
    uniquepos = SVector{3, T}[pge.pos[I[1]]]
    if vtypes !== nothing
        vtypeinit = (vtypes[I[1]], 1)
        symb_to_int = Dict{Vector{Tuple{valtype(vtypes),Int}}, Cint}([vtypeinit] => 1)
        last_types = [vtypeinit]
        types = Cint[1]
        new_symbol = 2
    else
        types = zeros(Cint, length(pge.pos))
    end
    for i in 2:length(pge.pos)
        j = I[i]
        if T <: Rational ? pge.pos[j] != last(uniquepos) : !isapprox(pge.pos[j], last(uniquepos); atol=0.0001)
            push!(uniquepos, pge.pos[j])
            if vtypes !== nothing
                last_types = [(vtypes[j], 1)]
                symb = get!(symb_to_int, last_types, new_symbol)
                push!(types, symb)
                new_symbol += (symb == new_symbol)
            end
        elseif vtypes !== nothing
            flag = true
            vtypej = vtypes[j]
            for k in 1:length(last_types)
                thistype, x = last_types[k]
                if thistype == vtypej
                    last_types[k] = (thistype, x+1)
                    flag = false
                    break
                end
            end
            if flag
                push!(last_types, (vtypej, 1))
            end
            symb = get!(symb_to_int, sort!(last_types), new_symbol)
            types[end] = symb
            new_symbol += (symb == new_symbol)
        end
    end

    n = length(uniquepos)
    positions = Matrix{Cdouble}(undef, 3, n)
    den = 1
    for i in 1:n
        if T <: Rational
            den = lcm(den, lcm(denominator.(uniquepos[i])))
        end
        positions[:,i] .= uniquepos[i]
    end

    ϵ = T <: Rational ? 100*eps(Cdouble) : 0.0001
    rotations = Array{Cint}(undef, 3, 3, 384)
    translations = Array{Cdouble}(undef, 3, 384)
    len = ccall((:spg_get_symmetry, libsymspg), Cint,
            (Ptr{Cint}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
            rotations, translations, 384, lattice, positions, types, n, ϵ)
    if T <: Rational
        den = lcm(len, den)
    end
    # The first symmetry should always be the identity, which we can skip
    rots = SMatrix{3,3,Int,9}[]
    transs = SVector{3,T}[]
    vmaps = Vector{PeriodicVertex3D}[]
    floatpos = [float(x) for x in pge.pos]
    hasmirror = false # whether a mirror plane exists or not. If so, the graph is achiral
    for i in 2:len
        rot = SMatrix{3,3,Cint,9}(transpose(rotations[:,:,i]))
        tr = SVector{3,Cdouble}(translations[:,i])
        trans = T <: Rational ? SVector{3,T}(round.(__rattype(T), den .* tr) .// den) : SVector{3,T}(tr)
        vmap = check_symmetry(pge, trans, rot, vtypes)
        if isnothing(vmap)
            trans = SVector{3,T}(pge.pos[last(findmin([norm(x .- tr) for x in floatpos]))])
            vmap = check_symmetry(pge, trans, rot, vtypes)
            isnothing(vmap) && continue
        end
        # if rot == LinearAlgebra.I
        #     error("The periodic graph is not minimal") # TODO: update pge?
        # end
        hasmirror |= det(rot) < 0
        push!(rots, rot)
        push!(vmaps, vmap::Vector{PeriodicVertex3D})
        push!(transs, trans)
    end
    return SymmetryGroup3D{T}(vmaps, rots, transs, hasmirror, length(pge.pos))
end
