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
