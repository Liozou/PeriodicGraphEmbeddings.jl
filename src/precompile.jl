macro enforce(expr) # strong @assert
    msg = string(expr)
    return :($(esc(expr)) ? $(nothing) : throw(AssertionError($msg)))
end

function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing

    # utils.jl
    @enforce precompile(Tuple{typeof(prepare_periodic_distance_computations), Matrix{Float64}})
    @enforce precompile(Tuple{typeof(prepare_periodic_distance_computations), SMatrix{3,3,Float64,9}})
    @enforce precompile(Tuple{typeof(periodic_distance!), MVector{3,Float64}, SVector{3,Float64}, Matrix{Float64}, Bool, Float64})
    @enforce precompile(Tuple{typeof(periodic_distance!), MVector{3,Float64}, SVector{3,Float64}, SMatrix{3,3,Float64,9}, Bool, Float64})
    @enforce precompile(Tuple{typeof(periodic_distance!), MVector{3,Float64}, MVector{3,Float64}, Matrix{Float64}, Bool, Float64})
    @enforce precompile(Tuple{typeof(periodic_distance!), MVector{3,Float64}, MVector{3,Float64}, SMatrix{3,3,Float64,9}, Bool, Float64})
    @enforce precompile(Tuple{typeof(periodic_distance!), MVector{3,Float64}, Matrix{Float64}, Bool, Float64})
    @enforce precompile(Tuple{typeof(periodic_distance!), MVector{3,Float64}, SMatrix{3,3,Float64,9}, Bool, Float64})
    @enforce precompile(Tuple{typeof(periodic_distance), SVector{3,Float64}, Matrix{Float64}, Bool, Float64})
    @enforce precompile(Tuple{typeof(periodic_distance), SVector{3,Float64}, Matrix{Float64}, Nothing, Nothing})
    @enforce precompile(Tuple{typeof(periodic_distance), SVector{3,Float64}, Matrix{Float64}})
    @enforce precompile(Tuple{typeof(periodic_distance), SVector{3,Float64}, SMatrix{3,3,Float64,9}, Bool, Float64})
    @enforce precompile(Tuple{typeof(periodic_distance), SVector{3,Float64}, SMatrix{3,3,Float64,9}, Nothing, Nothing})
    @enforce precompile(Tuple{typeof(periodic_distance), SVector{3,Float64}, SMatrix{3,3,Float64,9}})

    # types.jl
    @enforce precompile(Tuple{Type{EquivalentPosition}, SMatrix{3,3,Int,9}, SVector{3,Rational{Int}}})
    @enforce precompile(Tuple{Type{EquivalentPosition}, SMatrix{3,3,Int,9}, SVector{3,Float64}})
    @enforce precompile(Tuple{typeof(find_refid), Vector{String}})
    @enforce precompile(Tuple{typeof(cell_parameters), SMatrix{3,3,BigFloat,9}})
    for T in (Float64,Float32,Rational{Int})
        @enforce precompile(Tuple{Type{EquivalentPosition{T}}})
        @enforce precompile(Tuple{EquivalentPosition{T}, SMatrix{3,3,Int,9}})
        @enforce precompile(Tuple{EquivalentPosition{T}, SVector{3,Float64}})
        @enforce precompile(Tuple{EquivalentPosition{T}, SVector{3,Rational{Int}}})
        @enforce precompile(Tuple{typeof(Base.deepcopy_internal),Vector{EquivalentPosition{T}},IdDict{Any, Any}})
    end
    @enforce precompile(Tuple{typeof(parse), Type{EquivalentPosition}, String, NTuple{3,String}})
    @enforce precompile(Tuple{typeof(parse), Type{EquivalentPosition}, String})
    for T in (Float64,Float32,Rational{Int})
        @enforce precompile(Tuple{typeof(PeriodicGraphEmbeddings.__tostring), T, Bool, Bool})
        @enforce precompile(Tuple{typeof(show), Base.TTY, EquivalentPosition{T}})
        @enforce precompile(Tuple{typeof(show), Base.IOStream, EquivalentPosition{T}})
        @enforce precompile(Tuple{Type{Cell{T}}, Int, Matrix{BigFloat}, Vector{EquivalentPosition{T}}})
        @enforce precompile(Tuple{Type{Cell}, Int, Matrix{BigFloat}, Vector{EquivalentPosition{T}}})
        @enforce precompile(Tuple{Type{Cell}, Int, SMatrix{3,3,T,9}, Vector{EquivalentPosition{T}}})
        @enforce precompile(Tuple{Type{Cell}, Int, SMatrix{3,3,BigFloat,9}, Vector{EquivalentPosition{T}}})
        @enforce precompile(Tuple{typeof(cell_parameters), Cell{T}})
        @enforce precompile(Tuple{typeof(==), Cell{T}, Cell{T}})
        @enforce precompile(Tuple{Type{Cell{T}}, Int, NTuple{3,BigFloat}, NTuple{3,BigFloat}, Vector{EquivalentPosition{T}}})
        @enforce precompile(Tuple{Type{Cell{T}}, Int, NTuple{3,Int}, NTuple{3,Int}, Vector{EquivalentPosition{T}}})
        @enforce precompile(Tuple{Type{Cell}, Int, NTuple{3,BigFloat}, NTuple{3,BigFloat}, Vector{EquivalentPosition{T}}})
        @enforce precompile(Tuple{Type{Cell}, Int, NTuple{3,Int}, NTuple{3,Int}, Vector{EquivalentPosition{T}}})
        @enforce precompile(Tuple{Type{Cell{T}}, Int, NTuple{3,BigFloat}, NTuple{3,BigFloat}})
        @enforce precompile(Tuple{Type{Cell{T}}, Int, NTuple{3,Int}, NTuple{3,Int}})
        @enforce precompile(Tuple{Type{Cell{T}}})
        @enforce precompile(Tuple{Type{Cell}, Cell{T}, SMatrix{3,3,BigFloat,9}})
        @enforce precompile(Tuple{Type{Cell{T}}, SMatrix{3,3,BigFloat,9}})
        @enforce precompile(Tuple{typeof(show), Base.TTY, Cell{T}})
        @enforce precompile(Tuple{typeof(show), Base.IOStream, Cell{T}})
        @enforce precompile(Tuple{typeof(show), Base.TTY, MIME"text/plain", Cell{T}})
        @enforce precompile(Tuple{typeof(show), Base.IOStream, MIME"text/plain", Cell{T}})
    end
    @enforce precompile(Tuple{Type{Cell}, Int, NTuple{3,BigFloat}, NTuple{3,BigFloat}})
    @enforce precompile(Tuple{Type{Cell}})
    @enforce precompile(Tuple{Type{Cell}, SMatrix{3,3,BigFloat,9}})

    cell = Cell{Rational{Int}}
    for T in (Float32, Float64, Rational{Int32}, Rational{Int64}, Rational{Int128}, Rational{BigInt})
        for D in 0:3
            @enforce precompile(Tuple{Type{PeriodicGraphEmbedding{D,T}}, PeriodicGraph{D}, Vector{SVector{D,T}}, cell})
            @enforce precompile(Tuple{typeof(==), PeriodicGraphEmbedding{D,T}, PeriodicGraphEmbedding{D,T}})
            @enforce precompile(Tuple{Type{PeriodicGraphEmbedding{D,T}}, PeriodicGraph{D}, Vector{SVector{D,T}}})
            @enforce precompile(Tuple{Type{PeriodicGraphEmbedding{D}}, PeriodicGraph{D}, Vector{SVector{D,T}}, cell})
            @enforce precompile(Tuple{Type{PeriodicGraphEmbedding}, PeriodicGraph{D}, Vector{SVector{D,T}}, cell})
            @enforce precompile(Tuple{Type{PeriodicGraphEmbedding{D,T}}, cell})
            @enforce precompile(Tuple{typeof(PeriodicGraphEmbeddings.return_to_box), PeriodicGraph{D}, Matrix{T}})
            @enforce precompile(Tuple{Type{PeriodicGraphEmbedding{D,T}}, PeriodicGraph{D}, Matrix{T}, cell})
            @enforce precompile(Tuple{typeof(sort_periodicgraphembedding!), PeriodicGraph{D}, Matrix{T}, cell})
            @enforce precompile(Tuple{typeof(sort_periodicgraphembedding!), PeriodicGraph{D}, Matrix{T}})
            @enforce precompile(Tuple{typeof(getindex), PeriodicGraphEmbedding{D,T}, Int})
            @enforce precompile(Tuple{typeof(getindex), PeriodicGraphEmbedding{D,T}, PeriodicVertex{D}})
            @enforce precompile(Tuple{typeof(length), PeriodicGraphEmbedding{D,T}})
            @enforce precompile(Tuple{typeof(getindex), PeriodicGraphEmbedding{D,T}, Vector{Int}})
            for N in 0:3
                for S in (Float32, Float64, Rational{Int32}, Rational{Int64}, Rational{Int128}, Rational{BigInt})
                    @enforce precompile(Tuple{Type{PeriodicGraphEmbedding{D,T}}, PeriodicGraphEmbedding{N,S}})
                end
                @enforce precompile(Tuple{Type{PeriodicGraphEmbedding{D}}, PeriodicGraphEmbedding{N,T}})
            end
        end
    end

    for T in (Float64,Float32,Rational{Int})
        @enforce precompile(Tuple{typeof(getindex), Type{PeriodicSymmetry3D{T}}, Int})
        @enforce precompile(Tuple{typeof(getindex), Type{PeriodicSymmetry3D{T}}, PeriodicVertex3D})
        for S in (Float64,Float32,Rational{Int},Int)
            @enforce precompile(Tuple{PeriodicSymmetry3D{T}, Vector{S}})
            @enforce precompile(Tuple{PeriodicSymmetry3D{T}, SMatrix{3,3,S,9}})
            @enforce precompile(Tuple{PeriodicSymmetry3D{T}, Matrix{T}})
        end
        @enforce precompile(Tuple{SymmetryGroup3D{T}, Int})
        @enforce precompile(Tuple{typeof(unique), SymmetryGroup3D{T}})
        @enforce precompile(Tuple{typeof(getindex), SymmetryGroup3D{T}, Int})
        @enforce precompile(Tuple{typeof(iterate), SymmetryGroup3D{T}, Int})
        @enforce precompile(Tuple{typeof(iterate), SymmetryGroup3D{T}})
        @enforce precompile(Tuple{typeof(length), SymmetryGroup3D{T}})
        @enforce precompile(Tuple{typeof(one), SymmetryGroup3D{T}})
        @enforce precompile(Tuple{Type{SymmetryGroup3D}, Vector{Vector{PeriodicVertex3D}}, Vector{EquivalentPosition{T}}, Bool, Int})
        @enforce precompile(Tuple{Type{SymmetryGroup3D{T}}, Vector{Vector{PeriodicVertex3D}}, Vector{SMatrix{3,3,Int,9}}, Vector{SVector{3,T}}, Bool, Int})
    end

    # other.jl
    for D in 1:3
        for T in (Float64,Float32,Rational{Int})
            @enforce precompile(Tuple{typeof(PeriodicGraphs.offset_representatives!), PeriodicGraphEmbedding{D,T}, Vector{SVector{D,Int}}})
            @enforce precompile(Tuple{typeof(PeriodicGraphs.swap_axes!), PeriodicGraphEmbedding{D,T}, SVector{D,Int}})
            @enforce precompile(Tuple{typeof(PeriodicGraphs.swap_axes!), PeriodicGraphEmbedding{D,T}, NTuple{D,Int}})
            @enforce precompile(Tuple{typeof(PeriodicGraphs.swap_axes!), PeriodicGraphEmbedding{D,T}, Vector{Int}})
        end
    end

    # symmetries.jl
    @enforce precompile(Tuple{typeof(PeriodicGraphEmbeddings.joinletters), NTuple{11,Cchar}})
    @enforce precompile(Tuple{typeof(PeriodicGraphEmbeddings.joinletters), NTuple{17,Cchar}})
    @enforce precompile(Tuple{typeof(PeriodicGraphEmbeddings.joinletters), NTuple{6,Cchar}})
    @enforce precompile(Tuple{typeof(PeriodicGraphEmbeddings.joinletters), NTuple{7,Cchar}})
    @enforce precompile(Tuple{typeof(getproperty), PeriodicGraphEmbeddings.SpglibDataset, Symbol})
    @enforce precompile(Tuple{typeof(find_hall_number), String, String, Int, Bool})
    @enforce precompile(Tuple{typeof(find_hall_number), String, String, Int})
    @enforce precompile(Tuple{typeof(find_hall_number), String, String})
    @enforce precompile(Tuple{typeof(find_hall_number), String})
    @enforce precompile(Tuple{typeof(get_symmetry_equivalents), Rational{Int}, Int})
    @enforce precompile(Tuple{typeof(get_symmetry_equivalents), Float64, Int})
    @enforce precompile(Tuple{typeof(get_symmetry_equivalents), Float32, Int})
    @enforce precompile(Tuple{typeof(get_symmetry_equivalents), Int})
    for T in (Float32, Float64, Rational{Int32}, Rational{Int64}, Rational{Int128}, Rational{BigInt})
        @enforce precompile(Tuple{typeof(PeriodicGraphEmbeddings.get_spglib_dataset), PeriodicGraphEmbedding3D{T}, Nothing})
        @enforce precompile(Tuple{typeof(PeriodicGraphEmbeddings.get_spglib_dataset), PeriodicGraphEmbedding3D{T}})
        @enforce precompile(Tuple{typeof(PeriodicGraphEmbeddings.get_spglib_dataset), PeriodicGraphEmbedding3D{T}, Vector{Int}})
        @enforce precompile(Tuple{typeof(PeriodicGraphEmbeddings.get_spglib_dataset), PeriodicGraphEmbedding3D{T}, Vector{Symbol}})
        @enforce precompile(Tuple{typeof(check_valid_symmetry), PeriodicGraphEmbedding3D{T}, SVector{3,T}, Nothing, Nothing, Bool})
        @enforce precompile(Tuple{typeof(check_valid_symmetry), PeriodicGraphEmbedding3D{T}, SVector{3,T}, Nothing, Nothing})
        @enforce precompile(Tuple{typeof(check_valid_symmetry), PeriodicGraphEmbedding3D{T}, SVector{3,T}, Nothing})
        @enforce precompile(Tuple{typeof(check_valid_symmetry), PeriodicGraphEmbedding3D{T}, SVector{3,T}})
        @enforce precompile(Tuple{typeof(check_valid_symmetry), PeriodicGraphEmbedding3D{T}, SVector{3,T}, SMatrix{3,3,Int,9}})
        @enforce precompile(Tuple{typeof(check_valid_symmetry), PeriodicGraphEmbedding3D{T}, SVector{3,T}, SMatrix{3,3,Int,9}, Vector{Symbol}})
        @enforce precompile(Tuple{typeof(check_valid_symmetry), PeriodicGraphEmbedding3D{T}, SVector{3,T}, Nothing, Vector{Symbol}})
        @enforce precompile(Tuple{typeof(find_symmetries), PeriodicGraphEmbedding3D{T}, Nothing, typeof(check_valid_symmetry)})
        @enforce precompile(Tuple{typeof(find_symmetries), PeriodicGraphEmbedding3D{T}, Nothing})
        @enforce precompile(Tuple{typeof(find_symmetries), PeriodicGraphEmbedding3D{T}})
        @enforce precompile(Tuple{typeof(find_symmetries), PeriodicGraphEmbedding3D{T}, Vector{Symbol}, typeof(check_valid_symmetry)})
        @enforce precompile(Tuple{typeof(find_symmetries), PeriodicGraphEmbedding3D{T}, Vector{Symbol}})
    end

    # io.jl
    for T in (Float32, Float64, Rational{Int32}, Rational{Int64}, Rational{Int128}, Rational{BigInt})
        @enforce precompile(Tuple{typeof(export_vtf), String, PeriodicGraphEmbedding3D{T}, Nothing, Int, Bool})
        @enforce precompile(Tuple{typeof(export_vtf), String, PeriodicGraphEmbedding3D{T}, Nothing, Int})
        @enforce precompile(Tuple{typeof(export_vtf), String, PeriodicGraphEmbedding3D{T}, Nothing})
        @enforce precompile(Tuple{typeof(export_vtf), String, PeriodicGraphEmbedding3D{T}, Vector{Symbol}, Int, Bool})
        @enforce precompile(Tuple{typeof(export_vtf), String, PeriodicGraphEmbedding3D{T}, Vector{Symbol}, Int})
        @enforce precompile(Tuple{typeof(export_vtf), String, PeriodicGraphEmbedding3D{T}, Vector{Symbol}})
        @enforce precompile(Tuple{typeof(export_vtf), String, PeriodicGraphEmbedding3D{T}})
        for D in 0:3
            @enforce precompile(Tuple{typeof(export_cgd), String, PeriodicGraphEmbedding{D,T}, String, Bool})
            @enforce precompile(Tuple{typeof(export_cgd), String, PeriodicGraphEmbedding{D,T}, String})
            @enforce precompile(Tuple{typeof(export_cgd), String, PeriodicGraphEmbedding{D,T}})
        end
    end
    for D in 0:3
        @enforce precompile(Tuple{typeof(export_cgd), String, PeriodicGraph{D}, String, Bool})
        @enforce precompile(Tuple{typeof(export_cgd), String, PeriodicGraph{D}, String})
        @enforce precompile(Tuple{typeof(export_cgd), String, PeriodicGraph{D}})
    end
end

_precompile_()
