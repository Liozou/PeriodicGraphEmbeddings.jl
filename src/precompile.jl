using PeriodicGraphEmbeddings, PeriodicGraphs, Graphs, StaticArrays
using PrecompileTools

@static if VERSION < v"1.7-"
    Returns(a) = (_ -> a)
end

macro _precompile_pge(N, _name)
    g = Symbol(_name, :_g)
    v = Symbol(:onev, N)
    symb = gensym()
    ret = quote end
    for (T, _pos) in [(Float64, :_float),
                      (Rational{Int64}, :_rat64),
                      (Rational{Int128}, :_rat128),
                      (Rational{BigInt}, :_ratbig)
                     ]
        pos = Symbol(_name, :_pos, _pos)
        stacked = Symbol(pos, :_stacked)
        name = Symbol(_name, _pos)
        append!(ret.args, (quote
            $name = PeriodicGraphEmbedding($g, $pos)
            export_cgd(t, $name)
            export_cgd(t, $g)
            $name[1]
            $name[$v]
            $name[[2,1]]
            make_supercell($name, $(ntuple(Returns(2), N)))
            make_supercell($name, $(SVector{N,Int}(2 for _ in 1:N)))
            swap_axes!($name, $([i for i in 1:N]))
            swap_axes!($name, $(SVector{N,Int}([i for i in 1:N])))
            offset_representatives!($name, [SVector{$N,Int}(i+j*(-1)^i for j in 1:$N) for i in 1:length($name)])
            PeriodicGraphEmbedding($g, $stacked) == PeriodicGraphEmbedding($g, $stacked, cell)
            PeriodicGraphEmbedding{$N}($g, $stacked)
            PeriodicGraphEmbedding{$N}($g, $pos)
            SortedPeriodicGraphEmbedding($g, $stacked)
            SortedPeriodicGraphEmbedding($g, $stacked, cell)
            SortedPeriodicGraphEmbedding{$T}($g, $stacked)
        end).args)
        for i in N:3
            push!(ret.args, :($symb = PeriodicGraphEmbedding{$i}($name)))
            for j in N:(i-1)
                push!(ret.args, :(PeriodicGraphEmbedding{$j}($symb)))
            end
        end
    end
    return esc(ret)
end

@setup_workload begin
    apd_g = PeriodicGraph("3 1 3 0 0 0 1 44 0 0 0 2 5 1 0 0 2 21 0 0 0 2 84 0 0 0 2 93 0 0 0 3 46 0 0 0 3 50 0 0 0 3 79 0 0 0 4 6 0 0 0 4 43 0 0 0 4 60 0 0 0 4 74 0 0 0 5 78 0 0 0 6 17 0 0 0 7 69 0 0 0 7 73 0 0 0 8 35 0 0 0 8 52 0 0 0 8 62 0 0 0 8 76 0 0 0 9 16 0 0 0 9 22 0 0 0 9 30 1 0 0 9 86 0 0 0 10 13 0 0 0 10 54 0 0 0 10 85 0 0 0 10 87 0 1 0 11 45 1 0 0 11 58 0 0 0 11 81 0 0 0 11 92 0 0 0 12 40 0 0 0 12 69 0 0 0 13 77 0 0 0 14 26 0 0 -1 14 66 0 0 0 15 51 0 0 0 15 61 0 0 0 16 55 0 0 -1 17 31 0 0 0 17 35 0 0 0 17 96 0 0 0 18 24 0 0 0 18 28 0 0 0 18 50 0 0 0 18 72 0 0 1 19 51 0 0 0 19 55 0 0 0 20 36 1 0 0 20 68 0 0 1 20 82 0 0 0 20 95 0 0 0 21 26 0 0 0 22 66 0 0 0 23 48 0 0 0 23 83 0 0 0 23 93 0 0 0 23 96 1 0 0 24 57 -1 0 0 25 29 0 0 0 25 38 0 0 0 26 71 0 0 0 26 91 0 0 0 27 29 0 0 0 27 37 0 0 0 28 70 0 0 0 29 49 1 0 0 29 94 0 0 0 30 53 0 0 0 31 78 0 0 0 32 41 0 0 0 32 42 0 0 0 32 64 0 0 0 32 72 0 0 0 33 37 0 0 0 33 57 0 0 -1 34 61 0 0 0 34 69 0 0 0 36 73 0 0 0 37 42 0 0 0 37 56 0 0 0 38 65 0 1 0 38 67 0 0 0 38 85 0 0 0 39 44 0 0 0 39 57 0 0 0 40 65 0 0 0 40 81 0 0 0 40 95 0 0 0 41 77 0 0 0 43 53 0 0 0 44 84 0 0 0 44 94 0 0 0 45 61 0 0 0 46 77 0 0 0 47 51 0 0 0 47 73 0 0 0 48 59 0 0 0 49 77 0 0 0 51 75 0 0 1 52 61 0 0 0 53 75 0 0 0 53 76 0 0 0 54 88 0 0 0 55 82 0 0 0 55 92 0 0 0 56 89 0 0 0 57 91 0 0 0 58 59 0 0 0 59 62 0 0 0 59 86 0 0 0 60 70 0 0 -1 63 70 0 0 0 63 78 0 0 0 64 88 0 0 0 66 74 0 0 0 66 83 0 0 0 67 89 0 0 0 68 89 0 -1 0 69 87 0 0 0 70 71 0 0 0 73 90 0 0 1 78 79 0 0 0 80 88 0 0 0 80 89 -1 0 0 88 90 0 1 0")
    apd_pos_float = SVector{3, Float64}[[0.5, 0.7027, 0.6094], [0.8185, 0.5552, 0.6114], [0.316, 0.6986, 0.6073], [0.3185, 0.4448, 0.1114], [0.0, 0.5673, 0.6258], [0.25, 0.4698, 0.25], [0.25, 0.0302, 0.75], [0.316, 0.3014, 0.3927], [0.816, 0.3014, 0.1073], [0.3185, 0.9448, 0.3886], [0.816, 0.1986, 0.6073], [0.5, 0.0673, 0.6258], [0.2366, 0.876, 0.4299], [0.7137, 0.5, 0.0], [0.25, 0.218, 0.75], [0.75, 0.25, 0.0], [0.1815, 0.4448, 0.3886], [0.184, 0.6986, 0.8927], [0.5, 0.2027, 0.8906], [0.8185, 0.0552, 0.8886], [0.75, 0.5302, 0.75], [0.7634, 0.376, 0.0701], [0.8185, 0.4448, 0.3886], [0.0, 0.7027, 0.8906], [0.7634, 0.876, 0.4299], [0.6815, 0.5552, 0.8886], [0.75, 0.782, 0.25], [0.2366, 0.624, 0.9299], [0.816, 0.8014, 0.3927], [0.0, 0.2973, 0.1094], [0.2137, 0.5, 0.5], [0.316, 0.8014, 0.1073], [0.75, 0.75, 0.0], [0.2366, 0.124, 0.5701], [0.2634, 0.376, 0.4299], [0.0, 0.0673, 0.8742], [0.684, 0.8014, 0.1073], [0.6815, 0.9448, 0.3886], [0.75, 0.718, 0.75], [0.6815, 0.0552, 0.6114], [0.25, 0.782, 0.25], [0.5, 0.7973, 0.1094], [0.2366, 0.376, 0.0701], [0.684, 0.6986, 0.6073], [0.0, 0.2027, 0.6094], [0.25, 0.75, 0.5], [0.2634, 0.124, 0.9299], [0.7366, 0.376, 0.4299], [0.0, 0.7973, 0.3906], [0.25, 0.718, 0.75], [0.316, 0.1986, 0.8927], [0.25, 0.25, 0.5], [0.184, 0.3014, 0.1073], [0.25, 0.9698, 0.25], [0.684, 0.1986, 0.8927], [0.7366, 0.876, 0.0701], [0.816, 0.6986, 0.8927], [0.75, 0.25, 0.5], [0.684, 0.3014, 0.3927], [0.2863, 0.5, 0.0], [0.184, 0.1986, 0.6073], [0.5, 0.2973, 0.3906], [0.25, 0.5302, 0.75], [0.2634, 0.876, 0.0701], [0.7137, 0.0, 0.5], [0.6815, 0.4448, 0.1114], [0.75, 0.9698, 0.25], [0.7863, 0.0, 0.0], [0.3185, 0.0552, 0.6114], [0.3185, 0.5552, 0.8886], [0.5, 0.5673, 0.8742], [0.25, 0.75, 0.0], [0.1815, 0.0552, 0.8886], [0.5, 0.4327, 0.1258], [0.25, 0.25, 0.0], [0.25, 0.282, 0.25], [0.184, 0.8014, 0.3927], [0.1815, 0.5552, 0.6114], [0.2634, 0.624, 0.5701], [0.0, 0.9327, 0.1258], [0.7634, 0.124, 0.5701], [0.7366, 0.124, 0.9299], [0.75, 0.4698, 0.25], [0.7366, 0.624, 0.5701], [0.5, 0.9327, 0.3742], [0.75, 0.282, 0.25], [0.2863, 0.0, 0.5], [0.1815, 0.9448, 0.1114], [0.8185, 0.9448, 0.1114], [0.2137, 0.0, 0.0], [0.7634, 0.624, 0.9299], [0.75, 0.218, 0.75], [0.7863, 0.5, 0.5], [0.75, 0.75, 0.5], [0.75, 0.0302, 0.75], [0.0, 0.4327, 0.3742]]
    apd_pos_float_stacked = reduce(hcat, apd_pos_float)
    apd_pos_rat64 = [rationalize.(Int64, x) for x in apd_pos_float]
    apd_pos_rat64_stacked = reduce(hcat, apd_pos_rat64)
    apd_pos_rat128 = SVector{3,Rational{Int128}}.(apd_pos_rat64)
    apd_pos_rat128_stacked = reduce(hcat, apd_pos_rat128)
    apd_pos_ratbig = SVector{3,Rational{BigInt}}.(apd_pos_rat64)
    apd_pos_ratbig_stacked = reduce(hcat, apd_pos_ratbig)
    onepos3 = SVector{3, Float64}(1,2,3)
    onev3 = PeriodicVertex(3, (1,0,-1))
    hcb_g = PeriodicGraph("2 1 2 -1 0 1 2 0 0 1 2 0 1")
    hcb_pos_rat64 = SVector{2,Rational{Int64}}[(0,0), (1//3,-1//3)]
    hcb_pos_rat64_stacked = reduce(hcat, hcb_pos_rat64)
    hcb_pos_float = SVector{2,Float64}.(hcb_pos_rat64)
    hcb_pos_float_stacked = reduce(hcat, hcb_pos_float)
    hcb_pos_rat128 = SVector{2,Rational{Int128}}.(hcb_pos_rat64)
    hcb_pos_rat128_stacked = reduce(hcat, hcb_pos_rat128)
    hcb_pos_ratbig = SVector{2,Rational{BigInt}}.(hcb_pos_rat64)
    hcb_pos_ratbig_stacked = reduce(hcat, hcb_pos_ratbig)
    onepos2 = SVector{2,Float64}(1,2)
    onev2 = PeriodicVertex(2, (1,0))
    t = tempname(; cleanup=false)
    types = [:Si for _ in 1:length(apd_pos_float)]
    @compile_workload begin
        string(parse(EquivalentPosition, " - y; x+0.3, z-x"))
        parse(EquivalentPosition, "a+c, - a +x, c", ("x","a","c"))(onepos3)
        PeriodicGraphEmbeddings.find_refid(["X+1/2,-A,Z","+X,+Z,+A"])
        cell = Cell(144, (7.0, 7.0, 7.0), (90.0, 90.0, 90.0))
        cell == Cell{Float64}(15, (4.0, 4.0, 4.0), (90.0, 90.0, 90.0))
        cell_parameters(cell.mat)
        find_hall_number("i -4 -2")

        for (T, pos) in [(Float64, apd_pos_float),
                         (Rational{Int64}, apd_pos_rat64),
                         (Rational{Int128}, apd_pos_rat128),
                         (Rational{BigInt}, apd_pos_ratbig)]
            apd = PeriodicGraphEmbedding3D(apd_g, pos)
            apd_symm = find_symmetries(apd)
            find_symmetries(apd, types)
            find_symmetries(apd, zeros(Int, length(apd)))
            export_vtf(t, apd)
            export_vtf(t, apd, types)
            apd_symm(6)
            apd_symm[3](2)
            unique(apd_symm)
            length(apd_symm)
            one(apd_symm)
            symm = collect(apd_symm)[2]
            symm(3)
            symm(onev3)
            symm([3//4, 2//3, 1//1])
            symm(SVector{3,T}(T[1//2, 6//5, -1//10]))
            symm(T[1 0 0;0 1 0;0 0 1])
            symm(SMatrix{3,3,T,9}(T[1 0 0;0 1 0;0 0 1]))
            PeriodicGraphEmbeddings.get_spglib_dataset(apd).rotations
            PeriodicGraphEmbeddings.get_spglib_dataset(apd, types).hall_symbol
            RingAttributions(apd.g, apd_symm)
            RingAttributions(apd.g, true, apd_symm)
            RingAttributions(apd.g, 6, apd_symm)
            RingAttributions(apd.g, true, 6, apd_symm)
        end

        @_precompile_pge 3 apd
        @_precompile_pge 2 hcb
    end
    rm(t*".cgd")
    rm(t*".vtf")
end

