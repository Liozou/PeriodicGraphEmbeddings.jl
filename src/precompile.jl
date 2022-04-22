# types.jl
precompile(Tuple{Type{CrystalNets.EquivalentPosition}, SMatrix{3,3,Int,9}, SVector{3,Rational{Int}}})
precompile(Tuple{typeof(CrystalNets.find_refid), Vector{String}})
precompile(Tuple{typeof(parse), Type{CrystalNets.EquivalentPosition}, String, NTuple{3,String}})
precompile(Tuple{typeof(CrystalNets.rationaltostring), Int, Bool, Bool})
precompile(Tuple{typeof(show), Base.TTY, CrystalNets.EquivalentPosition})
precompile(Tuple{typeof(show), Base.IOStream, CrystalNets.EquivalentPosition})
precompile(Tuple{Type{cell}, Int, mat, Vector{CrystalNets.EquivalentPosition}})
precompile(Tuple{Type{cell}, Int, Matrix{BigFloat}, Vector{CrystalNets.EquivalentPosition}})
precompile(Tuple{typeof(==), cell, cell})
precompile(Tuple{Type{cell}, Int, NTuple{3,BigFloat}, NTuple{3,BigFloat}, Vector{CrystalNets.EquivalentPosition}})
precompile(Tuple{Type{cell}, Int, NTuple{3,Int}, NTuple{3,Int}, Vector{CrystalNets.EquivalentPosition}})
precompile(Tuple{Type{cell}})
precompile(Tuple{typeof(CrystalNets.cell_parameters), mat})
precompile(Tuple{typeof(CrystalNets.cell_parameters), fmat})
precompile(Tuple{typeof(CrystalNets.cell_parameters), cell})
precompile(Tuple{Type{cell}, mat})
precompile(Tuple{Type{cell}, fmat})
precompile(Tuple{typeof(show), Base.TTY, cell})
precompile(Tuple{typeof(show), Base.IOStream, cell})
precompile(Tuple{typeof(show), Base.TTY, MIME"text/plain", cell})
precompile(Tuple{typeof(show), Base.IOStream, MIME"text/plain", cell})
precompile(Tuple{typeof(CrystalNets.prepare_periodic_distance_computations), mat})
precompile(Tuple{typeof(CrystalNets.prepare_periodic_distance_computations), fmat})
precompile(Tuple{typeof(CrystalNets.periodic_distance!), MVector{3,Float64}, _pos, fmat, Bool, Float64})


# symmetries.jl
precompile(Tuple{typeof(CrystalNets.joinletters), NTuple{11,Cchar}})
precompile(Tuple{typeof(CrystalNets.joinletters), NTuple{17,Cchar}})
precompile(Tuple{typeof(CrystalNets.joinletters), NTuple{6,Cchar}})
precompile(Tuple{typeof(CrystalNets.joinletters), NTuple{7,Cchar}})
precompile(Tuple{typeof(getproperty), CrystalNets.SpglibDataset, Symbol})
precompile(Tuple{typeof(CrystalNets.find_hall_number), String, String, String})
precompile(Tuple{typeof(CrystalNets.get_symmetry_equivalents), Int})
for T in inttypes
    precompile(Tuple{typeof(CrystalNets.get_spglib_dataset), cnet{3,T}})
    precompile(Tuple{typeof(CrystalNets.find_symmetries), cnet{3,T}, collisions})
end

# topology.jl (TOCHANGE)
precompile(Tuple{typeof(CrystalNets.check_valid_symmetry), cnet{D,T}, pos{D,T}, Vector{collision}, Nothing})
precompile(Tuple{typeof(CrystalNets.check_valid_symmetry), cnet{D,T}, pos{D,T}, Vector{collision}})
