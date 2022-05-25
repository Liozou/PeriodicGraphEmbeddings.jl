macro enforce(expr) # strong @assert
    msg = string(expr)
    return :($(esc(expr)) ? $(nothing) : throw(AssertionError($msg)))
end

function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing

    # types.jl
    @enforce precompile(Tuple{Type{EquivalentPosition}, SMatrix{3,3,Int,9}, SVector{3,Rational{Int}}})
    @enforce precompile(Tuple{Type{EquivalentPosition}, SMatrix{3,3,Int,9}, SVector{3,Float64}})
    @enforce precompile(Tuple{typeof(Base.deepcopy_internal),Vector{EquivalentPosition{T}},IdDict{Any, Any}})
    @enforce precompile(Tuple{typeof(find_refid), Vector{String}})
    @enforce precompile(Tuple{typeof(parse), Type{EquivalentPosition{T}}, String, NTuple{3,String}})
    @enforce precompile(Tuple{typeof(PeriodicGraphEmbeddings.rationaltostring), Int, Bool, Bool})
    @enforce precompile(Tuple{typeof(show), Base.TTY, EquivalentPosition{T}})
    @enforce precompile(Tuple{typeof(show), Base.IOStream, EquivalentPosition{T}})
    @enforce precompile(Tuple{Type{cell}, Int, mat, Vector{EquivalentPosition{Rational{Int}}}})
    @enforce precompile(Tuple{Type{cell}, Int, Matrix{BigFloat}, Vector{EquivalentPosition{Rational{Int}}}})
    @enforce precompile(Tuple{typeof(==), cell, cell})
    @enforce precompile(Tuple{Type{cell}, Int, NTuple{3,BigFloat}, NTuple{3,BigFloat}, Vector{EquivalentPosition{Rational{Int}}}})
    @enforce precompile(Tuple{Type{cell}, Int, NTuple{3,Int}, NTuple{3,Int}, Vector{EquivalentPosition{Rational{Int}}}})
    @enforce precompile(Tuple{Type{cell}})
    @enforce precompile(Tuple{typeof(cell_parameters), mat})
    @enforce precompile(Tuple{typeof(cell_parameters), fmat})
    @enforce precompile(Tuple{typeof(cell_parameters), cell})
    @enforce precompile(Tuple{Type{cell}, mat})
    @enforce precompile(Tuple{Type{cell}, fmat})
    @enforce precompile(Tuple{typeof(show), Base.TTY, cell})
    @enforce precompile(Tuple{typeof(show), Base.IOStream, cell})
    @enforce precompile(Tuple{typeof(show), Base.TTY, MIME"text/plain", cell})
    @enforce precompile(Tuple{typeof(show), Base.IOStream, MIME"text/plain", cell})
    @enforce precompile(Tuple{typeof(prepare_periodic_distance_computations), mat})
    @enforce precompile(Tuple{typeof(prepare_periodic_distance_computations), fmat})
    @enforce precompile(Tuple{typeof(periodic_distance!), MVector{3,Float64}, _pos, fmat, Bool, Float64})
    @enforce precompile(Tuple{typeof(periodic_distance!), MVector{3,Float64}, fmat, Bool, Float64})


    # symmetries.jl
    @enforce precompile(Tuple{typeof(PeriodicGraphEmbeddings.joinletters), NTuple{11,Cchar}})
    @enforce precompile(Tuple{typeof(PeriodicGraphEmbeddings.joinletters), NTuple{17,Cchar}})
    @enforce precompile(Tuple{typeof(PeriodicGraphEmbeddings.joinletters), NTuple{6,Cchar}})
    @enforce precompile(Tuple{typeof(PeriodicGraphEmbeddings.joinletters), NTuple{7,Cchar}})
    @enforce precompile(Tuple{typeof(getproperty), PeriodicGraphEmbeddings.SpglibDataset, Symbol})
    @enforce precompile(Tuple{typeof(find_hall_number), String, String, String})
    @enforce precompile(Tuple{typeof(get_symmetry_equivalents), Rational{Int}, Int})
    @enforce precompile(Tuple{typeof(get_symmetry_equivalents), Float64, Int})
    @enforce precompile(Tuple{typeof(get_symmetry_equivalents), Int})
    for T in inttypes
        @enforce precompile(Tuple{typeof(PeriodicGraphEmbeddings.get_spglib_dataset), cnet{3,T}})
        @enforce precompile(Tuple{typeof(find_symmetries), cnet{3,T}, collisions})
    end

    # topology.jl (TOCHANGE)
    @enforce precompile(Tuple{typeof(check_valid_symmetry), cnet{D,T}, pos{D,T}, Vector{collision}, Nothing})
    @enforce precompile(Tuple{typeof(check_valid_symmetry), cnet{D,T}, pos{D,T}, Vector{collision}})


    # io (TOCHANGE)
    @enforce @enforce precompile(Tuple{typeof(export_vtf), String, crystclust})
    @enforce @enforce precompile(Tuple{typeof(export_cgd), String, crystclust})
end

#_precompile_()
