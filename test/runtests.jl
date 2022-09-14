using Test
using PeriodicGraphEmbeddings, PeriodicGraphs, Graphs, StaticArrays

using Aqua
Aqua.test_all(PeriodicGraphs)

@testset "EquivalentPosition" begin
    @test PeriodicGraphEmbeddings.find_refid(String[]) == ("x", "y", "z")
    @test PeriodicGraphEmbeddings.find_refid(["X+1/2,-A,Z","+X,+Z,+A"]) == ("x", "z", "a")
    @test_throws ErrorException PeriodicGraphEmbeddings.find_refid(["a,b"])
    @test_throws ErrorException PeriodicGraphEmbeddings.find_refid(["x,y,z,t"])
    @test string(parse(EquivalentPosition, " - y; x+0.3, z-x")) == "-y,x+3/10,-x+z"
    @test string(parse(EquivalentPosition, "12/6, x+y- 1z ,1-0z")) == "2,x+y-z,1"
    @test string(parse(EquivalentPosition, "a+c, - a +x, c", ("x","a","c"))) == "y+z,x-y,z"
    @test string(parse(EquivalentPosition, "y-2y,,- x+0.5+1/3")) == "-y,0,-x+5/6"
    @test_throws ErrorException parse(EquivalentPosition, "x")
    @test_throws ErrorException parse(EquivalentPosition, "x,y,z,t")
    @test parse(EquivalentPosition, ",,x")([1,2,3]) == [0,0,1]
    @test parse(EquivalentPosition, "z,y+1,-x-0.5")([1,2,3]) == [3,3,-3//2]
end

@testset "Cell" begin
    @test Cell() == Cell{Rational{Int}}(1, one(SMatrix{3,3,BigFloat,9}).*10, EquivalentPosition{Rational{Int}}[])
    @test Cell{Float64}() == Cell{Float64}(1, (10, 10, 10), (90, 90, 90))
    cell = Cell(2, (3,4,5), (60,60,90))
    @test only(cell.equivalents)([6,8//7,9]) == [-6,-8//7,-9]
    ((a, b, c), (α, β, γ)), mat = cell_parameters(cell)
    @test all(isapprox.(Float32.((a, b, c, α, β, γ)), (3, 4, 5, 60, 60, 90)))
    @test isapprox(mat, cell.mat)
    shiftright = Int[0 0 1; 1 0 0; 0 1 0]
    (a, b, c), (α, β, γ) = cell_parameters(shiftright * cell.mat)
    @test all(isapprox.(Float32.((a, b, c, α, β, γ)), (3, 4, 5, 60, 60, 90)))
    @test string(cell) == "Cell(2, (3.0, 4.0, 5.0), (60.0, 60.0, 90.0))"
    @test occursin("-P 1 (triclinic)", repr(MIME"text/plain"(), cell))
    @test Cell(SMatrix{3,3,BigFloat,9}(shiftright)) == Cell{Rational{Int}}(one(SMatrix{3,3,BigFloat,9}))
end

@testset "PeriodicGraphEmbedding" begin
    g = PeriodicGraph2D("2  1 2 0 0  1 3 0 0  2 3 0 0  2 1 1 0  3 1 0 1  3 2 0 1")
    pos = SVector{2,Rational{Int}}[(1//4, 1//3), (3//4, 1//3), (1//2, 2//3)]
    pge = PeriodicGraphEmbedding(g, pos)
    @test pge == PeriodicGraphEmbedding{2}(g, pos) == PeriodicGraphEmbedding{2,Rational{UInt}}(g, pos)
    pos_t = [[1//4, 1//3], [3//4, 1//3], [1//2, 2//3]]
    @test pge == PeriodicGraphEmbedding(g, pos_t) ==  PeriodicGraphEmbedding{2,Rational{Int}}(g, pos_t)

    struct BadIterator; x; end
    Base.IteratorEltype(::Type{BadIterator}) = Base.EltypeUnknown()
    Base.IteratorSize(::Type{BadIterator}) = Base.SizeUnknown()
    Base.iterate(b::BadIterator) = iterate(b.x)
    Base.iterate(b::BadIterator, foo) = iterate(b.x, foo)
    @test collect(BadIterator([3,6,7])) == [3,6,7]
    @test PeriodicGraphEmbedding{2,Float64}(g, BadIterator(pos)) == PeriodicGraphEmbedding{2,Float64}(g, pos_t)
    @test_throws MethodError PeriodicGraphEmbedding(g, BadIterator(pos))

    ofs_g = PeriodicGraph2D("2  1 2 0 0  1 3 0 1  2 3 0 1  2 1 1 0  3 1 0 0  3 2 0 0")
    ofs_pos = [1//4 3//4 1//2; 1//3 1//3 -1//3]
    @test PeriodicGraphEmbedding(ofs_g, ofs_pos) == pge
    @test ofs_g == g

    @test !issorted([pge[1], pge[2], pge[3]])
    posmat = reduce(hcat, pos)
    new_pge, perm = SortedPeriodicGraphEmbedding{Rational{Int128}}(copy(g), posmat)
    @test new_pge isa PeriodicGraphEmbedding{2,Rational{Int128}}
    @test perm == [1, 3, 2]
    @test issorted([new_pge[1], new_pge[2], new_pge[3]])
    @test ofs_g == g
    @test new_pge == pge[perm]
    small_new_pge, small_new_perm = SortedPeriodicGraphEmbedding(copy(g), posmat)
    @test small_new_perm ==  perm
    @test small_new_pge == new_pge
    @test small_new_pge isa PeriodicGraphEmbedding{2,Rational{Int32}}
    @test SortedPeriodicGraphEmbedding(PeriodicGraph3D(), fill(0, 3, 0))[1] == PeriodicGraphEmbedding3D{Rational{Int32}}(PeriodicGraph3D(), [])
    bigrational = 1//(1+big(typemax(UInt128)))
    bigpos = fill(bigrational, 1, 1)
    big_pge, = SortedPeriodicGraphEmbedding(PeriodicGraph1D(1), bigpos)
    @test big_pge isa PeriodicGraphEmbedding{1,Rational{BigInt}}
    @test_throws InexactError SortedPeriodicGraphEmbedding{Rational{Int}}(PeriodicGraph1D(1), bigpos)
    @test only(only(SortedPeriodicGraphEmbedding{Float64}(PeriodicGraph1D(1), bigpos)[1].pos)) ≈ bigrational

    @test pge[PeriodicVertex2D(3, (0,-1))] == [1//2, -1//3]
    @test PeriodicGraphEmbedding{2,Rational{Int}}(PeriodicGraphEmbedding{3,Rational{Int8}}(pge)) == pge
    @test_throws DimensionMismatch PeriodicGraphEmbedding{1,Float64}(pge)
end

@testset "PeriodicSymmetry3D" begin
    g = PeriodicGraph2D("2  1 2 0 0  1 3 0 0  2 3 0 0  2 1 1 0  3 1 0 1  3 2 0 1")
    pos = SVector{2,Rational{Int}}[(1//4, 1//3), (3//4, 1//3), (1//2, 2//3)]
    pge = PeriodicGraphEmbedding(g, pos)
    symm3D = find_symmetries(PeriodicGraphEmbedding3D(pge), [Symbol("") for _ in 1:length(pge)])
    @test symm3D(1) == symm3D(2) != symm3D(3)
    @test unique(symm3D) == [1, 3]
    allsymm3D = collect(symm3D)
    @test allsymm3D == PeriodicGraphEmbeddings.PeriodicSymmetry3D{Rational{Int}}[symm3D[1], symm3D[2], symm3D[3]]
    @test one(symm3D) isa PeriodicGraphEmbeddings.PeriodicSymmetry3D{Rational{Int}}
    @test one(symm3D) ∉ allsymm3D
    for symm in symm3D
        @test symm(3) == 3
        for i in 1:3
            @test symm(symm(i)) == i
            @test symm(symm(PeriodicVertex(i, (6,7,8)))) == PeriodicVertex(i, (6,7,8))
        end
        @test iszero(symm(zero(SVector{3,Float64})))
    end
    @test first(symm3D)([1//2, 1//3, 0]) == [-1//2, 1//3, 0]
    @test first(symm3D)([1 0 0; 2 1 0; 0 0 0]) == [-1 0 0; 2 1 0; 0 0 0]
    @test only(find_symmetries(PeriodicGraphEmbedding3D(pge), [1, 2, 3]))([7,8,9]) == [7,8,-9]
end

@testset "General symmetries" begin
    mog = PeriodicGraphEmbedding(
        PeriodicGraph3D("3 1 4 0 -1 -1 1 4 0 0 -1 1 5 0 0 0 1 11 -1 0 0 2 3 0 -1 0 2 3 0 0 0 2 5 0 0 0 2 11 -1 0 0 3 6 0 0 0 3 12 -1 0 0 4 6 0 0 0 4 12 -1 0 0 5 7 0 0 0 5 8 0 0 0 6 9 0 0 0 6 10 0 0 0 7 10 0 -1 -1 7 10 0 0 -1 7 11 0 0 0 8 9 0 -1 0 8 9 0 0 0 8 11 0 0 0 9 12 0 0 0 10 12 0 0 0"),
        SVector{3, Rational{Int32}}[[0//1, 0//1, 0//1], [0//1, 0//1, 1//3], [0//1, 1//2, 1//2], [0//1, 1//2, 5//6], [1//4, 0//1, 1//6], [1//4, 1//2, 2//3], [1//2, 0//1, 0//1], [1//2, 0//1, 1//3], [1//2, 1//2, 1//2], [1//2, 1//2, 5//6], [3//4, 0//1, 1//6], [3//4, 1//2, 2//3]],
        Cell(1, (8.7371, 4.8692, 10.7217), (90.0, 90.193, 90.0))
    )
    mogsymms = find_symmetries(mog, nothing, (a,b,c,d) -> check_valid_symmetry(a,b,c,d, true))
    @test mogsymms.vmaps == find_symmetries(mog).vmaps
    @test length(mogsymms) == 31
    @test one(mogsymms) ∉ mogsymms
    @test length(one(mogsymms).vmap) == length(first(mogsymms).vmap) == length(mog)

    @test find_hall_number("-C 2 2") == find_hall_number("Cmmm") == find_hall_number("C m m m")
    @test 1 != find_hall_number("C 2 2 -1bc", "ccca") != find_hall_number("C 2 2 -1bc", "cccb") != 1
    @test 1 != find_hall_number("a 2 2 -1ac", "a c a a") != find_hall_number("a b a a") != 1
    @test 1 != find_hall_number("b 2 2 -1bc") != find_hall_number("b 2 2 -1bc", "bbab") != 1
    @test PeriodicGraphEmbeddings.RAW_SYMMETRY_DATA[find_hall_number("", "", 127)][2] == 127
    @test 1 == redirect_stderr(devnull) do
                     find_hall_number("-C a 2 m", "abca", 0, true)
                end
    
    spgdataset = PeriodicGraphEmbeddings.get_spglib_dataset(mog, [mogsymms(i) for i in 1:length(mog)])
    @test spgdataset.hall_symbol == PeriodicGraphEmbeddings.get_spglib_dataset(mog).hall_symbol
end

@testset "Complex modifications" begin
    itn = PeriodicGraphEmbedding(
        PeriodicGraph3D("3 1 28 0 0 0 1 29 0 0 0 1 77 0 0 1 1 157 -1 0 0 2 29 0 0 0 2 30 0 0 0 2 31 0 0 0 2 33 0 0 0 3 31 0 0 0 3 32 0 0 0 3 40 0 0 1 3 127 -1 0 0 4 33 0 0 0 4 34 0 0 0 4 35 0 0 0 4 36 0 0 0 5 36 0 0 0 5 37 0 0 0 5 38 0 0 0 5 79 0 1 0 6 38 0 0 0 6 39 0 0 0 6 43 0 0 0 6 47 0 0 0 7 39 0 0 0 7 40 0 0 0 7 41 0 0 0 7 42 0 0 0 8 43 0 -1 0 8 44 0 0 0 8 45 0 0 0 8 46 0 0 0 9 47 0 0 0 9 48 0 0 0 9 49 0 0 0 9 50 0 0 0 10 35 0 0 0 10 50 0 0 0 10 51 0 0 0 10 52 0 0 0 11 52 0 0 0 11 53 0 0 0 11 54 0 0 0 11 138 0 0 1 12 54 0 0 0 12 55 0 0 0 12 56 0 0 0 12 57 0 0 0 13 49 0 0 0 13 57 0 0 0 13 58 0 0 0 13 62 0 0 0 14 58 0 0 0 14 59 0 0 0 14 60 0 0 0 14 67 0 1 0 15 51 0 0 0 15 60 0 0 0 15 61 0 0 0 15 150 0 0 1 16 62 0 0 0 16 63 0 0 0 16 64 0 1 0 16 146 0 0 0 17 65 0 0 0 17 66 0 0 0 17 67 0 0 0 17 68 0 0 0 18 55 0 0 0 18 68 0 0 0 18 69 0 0 0 18 70 0 0 0 19 53 0 0 0 19 70 0 0 0 19 71 0 0 0 19 141 0 0 1 20 61 0 -1 0 20 71 0 0 0 20 72 0 0 0 20 73 0 0 0 21 48 0 -1 0 21 66 0 0 0 21 73 0 0 0 21 74 0 0 0 22 46 0 0 0 22 74 0 0 0 22 75 0 0 0 22 78 0 0 0 23 42 0 0 0 23 75 0 0 0 23 76 0 0 0 23 77 0 0 0 24 78 0 0 0 24 79 0 0 0 24 80 0 0 0 24 112 -1 0 0 25 34 0 -1 0 25 72 0 0 0 25 80 0 0 0 25 81 0 0 0 26 32 0 -1 0 26 81 0 0 0 26 82 0 0 0 26 119 -1 0 0 27 28 0 0 0 27 44 0 0 1 27 82 0 0 0 27 123 -1 0 0 30 106 -1 0 0 37 108 -1 0 0 41 109 -1 0 0 45 85 -1 0 0 56 93 0 0 1 59 101 0 0 1 64 98 0 0 0 65 98 0 0 0 69 97 0 0 1 76 83 -1 0 0 83 110 0 0 0 83 111 0 0 0 83 158 0 0 -1 84 111 0 0 0 84 112 0 0 0 84 113 0 0 0 84 115 0 0 0 85 113 0 0 0 85 114 0 0 0 85 122 0 0 -1 86 115 0 0 0 86 116 0 0 0 86 117 0 0 0 86 118 0 0 0 87 118 0 0 0 87 119 0 0 0 87 120 0 0 0 87 160 0 -1 0 88 120 0 0 0 88 121 0 0 0 88 125 0 0 0 88 129 0 0 0 89 121 0 0 0 89 122 0 0 0 89 123 0 0 0 89 124 0 0 0 90 125 0 1 0 90 126 0 0 0 90 127 0 0 0 90 128 0 0 0 91 129 0 0 0 91 130 0 0 0 91 131 0 0 0 91 132 0 0 0 92 117 0 0 0 92 132 0 0 0 92 133 0 0 0 92 134 0 0 0 93 134 0 0 0 93 135 0 0 0 93 136 0 0 0 94 136 0 0 0 94 137 0 0 0 94 138 0 0 0 94 139 0 0 0 95 131 0 0 0 95 139 0 0 0 95 140 0 0 0 95 144 0 0 0 96 140 0 0 0 96 141 0 0 0 96 142 0 0 0 96 148 0 -1 0 97 133 0 0 0 97 142 0 0 0 97 143 0 0 0 98 144 0 0 0 98 145 0 0 0 99 146 0 0 0 99 147 0 0 0 99 148 0 0 0 99 149 0 0 0 100 137 0 0 0 100 149 0 0 0 100 150 0 0 0 100 151 0 0 0 101 135 0 0 0 101 151 0 0 0 101 152 0 0 0 102 143 0 1 0 102 152 0 0 0 102 153 0 0 0 102 154 0 0 0 103 130 0 1 0 103 147 0 0 0 103 154 0 0 0 103 155 0 0 0 104 128 0 0 0 104 155 0 0 0 104 156 0 0 0 104 159 0 0 0 105 124 0 0 0 105 156 0 0 0 105 157 0 0 0 105 158 0 0 0 106 159 0 0 0 106 160 0 0 0 106 161 0 0 0 107 116 0 1 0 107 153 0 0 0 107 161 0 0 0 107 162 0 0 0 108 114 0 1 0 108 162 0 0 0 108 163 0 0 0 109 110 0 0 0 109 126 0 0 -1 109 163 0 0 0"),
        SVector{3, Float64}[[0.05179, 0.48234, 0.91511], [0.07258, 0.70899, 0.79527], [0.05733, 0.85943, 0.91674], [0.1567, 0.81073, 0.70182], [0.0769, 0.8981, 0.42328], [0.16828, 0.85774, 0.34649], [0.10153, 0.72914, 0.17864], [0.09116, 0.1094, 0.1678], [0.31108, 0.85323, 0.56472], [0.30251, 0.71347, 0.78004], [0.39421, 0.47891, 0.91981], [0.49287, 0.50757, 0.83426], [0.45775, 0.75419, 0.66662], [0.49044, 0.89125, 0.8331], [0.39214, 0.86637, 0.90875], [0.47586, 0.88681, 0.47318], [0.45639, 0.13865, 0.65411], [0.48678, 0.26845, 0.82588], [0.38987, 0.24291, 0.90373], [0.29803, 0.10165, 0.77347], [0.31039, 0.10906, 0.56289], [0.16553, 0.2431, 0.33869], [0.09973, 0.48917, 0.16751], [0.07315, 0.14459, 0.41455], [0.15287, 0.05735, 0.69252], [0.06908, 0.09844, 0.7854], [0.04591, 0.24206, 0.90486], [0.07263, 0.35707, 0.92175], [0.0725, 0.57489, 0.83972], [0.00058, 0.77428, 0.6721], [0.08608, 0.75961, 0.89576], [0.08104, 0.97465, 0.87782], [0.13104, 0.72591, 0.77318], [0.17104, 0.92004, 0.74264], [0.22511, 0.74569, 0.73294], [0.09939, 0.8519, 0.55992], [0.00434, 0.87057, 0.33859], [0.13213, 0.83719, 0.41563], [0.16005, 0.75958, 0.29081], [0.08502, 0.82171, 0.05585], [0.03519, 0.73136, 0.17922], [0.12528, 0.60347, 0.18977], [0.13276, 0.97973, 0.23891], [0.07184, 0.14826, 0.03337], [0.02374, 0.11752, 0.16456], [0.13671, 0.19115, 0.23369], [0.24679, 0.85725, 0.43681], [0.31033, 0.97727, 0.57537], [0.37877, 0.81206, 0.57793], [0.30901, 0.76198, 0.66915], [0.34995, 0.76744, 0.87998], [0.32681, 0.57668, 0.83769], [0.37424, 0.35587, 0.9237], [0.437, 0.50174, 0.86299], [0.50488, 0.39341, 0.81252], [0.56233, 0.51857, 0.94629], [0.46746, 0.61737, 0.71605], [0.49311, 0.79797, 0.77808], [0.55586, 0.85508, 0.96913], [0.42348, 0.89731, 0.83151], [0.34179, 0.97775, 0.87295], [0.49145, 0.79108, 0.59376], [0.39606, 0.91647, 0.3709], [0.5, 0.0, 0.5], [0.48478, 0.16087, 0.57224], [0.37646, 0.14476, 0.57126], [0.49196, 0.01363, 0.75649], [0.47231, 0.23381, 0.71548], [0.54809, 0.1761, 0.95049], [0.41965, 0.272, 0.82398], [0.3216, 0.19659, 0.83529], [0.21912, 0.10632, 0.71656], [0.31008, 0.12923, 0.66896], [0.24474, 0.18938, 0.43561], [0.15436, 0.37958, 0.28127], [0.02837, 0.48757, 0.15061], [0.0914, 0.48486, 0.04877], [0.12655, 0.21322, 0.40253], [0.07155, 0.03412, 0.37835], [0.09499, 0.10387, 0.55069], [0.12619, 0.09926, 0.75981], [0.07463, 0.1931, 0.84397], [0.94821, 0.51766, 0.08489000000000002], [0.92742, 0.29101, 0.20472999999999997], [0.94267, 0.14056999999999997, 0.08326], [0.8432999999999999, 0.18927000000000005, 0.29818], [0.9231, 0.10189999999999999, 0.57672], [0.83172, 0.14226000000000005, 0.65351], [0.89847, 0.27086, 0.82136], [0.90884, 0.8906000000000001, 0.8322], [0.68892, 0.14676999999999996, 0.43528], [0.6974899999999999, 0.28652999999999995, 0.21996000000000004], [0.60579, 0.52109, 0.08018999999999998], [0.5071300000000001, 0.49243000000000003, 0.16574], [0.54225, 0.24580999999999997, 0.33338], [0.50956, 0.10875000000000001, 0.16690000000000005], [0.6078600000000001, 0.13363000000000003, 0.09125000000000005], [0.52414, 0.11319000000000001, 0.5268200000000001], [0.5436099999999999, 0.8613500000000001, 0.34589000000000003], [0.51322, 0.7315499999999999, 0.17412000000000005], [0.6101300000000001, 0.75709, 0.09626999999999997], [0.70197, 0.89835, 0.22653], [0.6896100000000001, 0.89094, 0.43711], [0.83447, 0.7569, 0.6613100000000001], [0.90027, 0.51083, 0.83249], [0.92685, 0.85541, 0.58545], [0.8471299999999999, 0.94265, 0.30748], [0.93092, 0.90156, 0.2146], [0.95409, 0.7579400000000001, 0.09514], [0.92737, 0.64293, 0.07825000000000004], [0.9275, 0.42511, 0.16027999999999998], [0.99942, 0.22572000000000003, 0.32789999999999997], [0.91392, 0.24039, 0.10424], [0.91896, 0.025349999999999984, 0.12217999999999996], [0.86896, 0.27408999999999994, 0.22682000000000002], [0.82896, 0.07996000000000003, 0.25736000000000003], [0.77489, 0.25431000000000004, 0.26705999999999996], [0.90061, 0.1481, 0.44008], [0.99566, 0.12943000000000005, 0.66141], [0.86787, 0.16281, 0.5843700000000001], [0.83995, 0.24041999999999997, 0.70919], [0.91498, 0.17828999999999995, 0.94415], [0.96481, 0.26864, 0.8207800000000001], [0.8747199999999999, 0.39653000000000005, 0.81023], [0.86724, 0.02027000000000001, 0.76109], [0.92816, 0.8517399999999999, 0.96663], [0.97626, 0.88248, 0.83544], [0.86329, 0.8088500000000001, 0.76631], [0.7532099999999999, 0.14275000000000004, 0.5631900000000001], [0.68967, 0.022730000000000028, 0.42462999999999995], [0.62123, 0.18794, 0.42206999999999995], [0.69099, 0.23802, 0.33085], [0.65005, 0.23256, 0.12002000000000002], [0.67319, 0.42332000000000003, 0.16230999999999995], [0.62576, 0.64413, 0.07630000000000003], [0.563, 0.49826000000000004, 0.13700999999999997], [0.49512, 0.60659, 0.18747999999999998], [0.43767, 0.48143, 0.053710000000000035], [0.53254, 0.38263, 0.28395000000000004], [0.5068900000000001, 0.20203000000000004, 0.22192], [0.44414, 0.14492000000000005, 0.030869999999999953], [0.5765199999999999, 0.10268999999999995, 0.16849000000000003], [0.65821, 0.022249999999999992, 0.12705], [0.5085500000000001, 0.20892, 0.40624000000000005], [0.6039399999999999, 0.08353, 0.6291], [0.51522, 0.8391299999999999, 0.42776000000000003], [0.62354, 0.85524, 0.42874], [0.50804, 0.98637, 0.24351], [0.52769, 0.76619, 0.28452], [0.45191000000000003, 0.8239, 0.049510000000000054], [0.5803499999999999, 0.728, 0.17601999999999995], [0.6784, 0.80341, 0.16471000000000002], [0.78088, 0.89368, 0.28344], [0.68992, 0.87077, 0.33104], [0.75526, 0.81062, 0.56439], [0.84564, 0.62042, 0.71873], [0.97163, 0.5124299999999999, 0.84939], [0.9086, 0.5151399999999999, 0.95123], [0.8734500000000001, 0.78678, 0.59747], [0.92845, 0.96588, 0.62165], [0.90501, 0.89613, 0.44931], [0.87381, 0.90074, 0.24019000000000001], [0.92537, 0.8069, 0.15603]],
        Cell(1, (24.2005, 12.5399, 14.2502), (72.599, 123.128, 90.0))
    )
    n = length(itn)
    itncopy = copy(itn)
    @test itncopy == itn !== itncopy

    itnoffset = copy(itncopy)
    offset_representatives!(itnoffset, [(-i,2,i+1) for i in 1:n])
    @test itn == itncopy != itnoffset
    for i in 1:n
        @test Float32.(itn[i]) == Float32.(itnoffset[PeriodicVertex(i, (-i,2,i+1))])
    end

    itnrotated = copy(itncopy)
    swap_axes!(itnrotated, [2,3,1])
    @test itn == itncopy != itnrotated
    for i in 1:n
        @test itn[i][[2,3,1]] == itnrotated[i]
    end

    itn312 = make_supercell(itn, (3, 1, 2))
    @test itn312 == make_supercell(itn, SVector{3,Int}(3, 1, 2)) == make_supercell(itn, [3, 1, 2])
    mat = itn.cell.mat
    mat312 = itn312.cell.mat
    for p in PeriodicGraphs.MetaClock([3, 1, 2])
        factor = n*(p[1] + 3*p[2] + 3*p[3])
        for i in 1:n
            @test Float32.(mat*itn[PeriodicVertex3D(i, p)]) == Float32.(mat312*itn312[factor+i])
        end
    end
end
