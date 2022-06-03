var documenterSearchIndex = {"docs":
[{"location":"types/#Types","page":"Types","title":"Types","text":"","category":"section"},{"location":"types/#Manual","page":"Types","title":"Manual","text":"","category":"section"},{"location":"types/","page":"Types","title":"Types","text":"PeriodicGraphEmbeddings provide the new type PeriodicGraphEmbedding which wraps:","category":"page"},{"location":"types/","page":"Types","title":"Types","text":"A PeriodicGraph\nThe list of positions of the vertices in a unit cell of the graph\nOptionally, a Cell if the dimension of the graph is 3 or below, which contains the geometry of the unit cell.","category":"page"},{"location":"types/","page":"Types","title":"Types","text":"A PeriodicGraphEmbedding can be built through different methods, depending on whether the list of positions should be permuted to be sorted, or offset to have all positions between 0 and 1 for instance:","category":"page"},{"location":"types/","page":"Types","title":"Types","text":"PeriodicGraphEmbedding\nPeriodicGraphEmbedding3D\nPeriodicGraphEmbedding{D,T}(graph::PeriodicGraph{D}, placement::AbstractMatrix{T}, cell::Cell) where {D,T}\nSortedPeriodicGraphEmbedding\nSortedPeriodicGraphEmbedding{T}(graph::PeriodicGraph{D}, placement::AbstractMatrix, cell::Cell) where {D,T}\nPeriodicGraphEmbedding{D,T}(pge::PeriodicGraphEmbedding{N,S}) where {D,T,N,S}","category":"page"},{"location":"types/#PeriodicGraphEmbeddings.PeriodicGraphEmbedding","page":"Types","title":"PeriodicGraphEmbeddings.PeriodicGraphEmbedding","text":"PeriodicGraphEmbedding{D,T}\n\nEmbedding in euclidean space of a PeriodicGraph of dimension D. Each vertex is assigned a D-uplet of coordinates of type T.\n\nPeriodicGraphEmbedding3D is provided as an alias for PeriodicGraphEmbedding{3}. Symmetry detection provided by PeriodicGraphEmbeddings.jl can only be performed on PeriodicGraphEmbedding3D.\n\n\n\n\n\n","category":"type"},{"location":"types/#PeriodicGraphEmbeddings.PeriodicGraphEmbedding3D","page":"Types","title":"PeriodicGraphEmbeddings.PeriodicGraphEmbedding3D","text":"PeriodicGraphEmbedding3D\n\nAlias for PeriodicGraphEmbedding{3}\n\n\n\n\n\n","category":"type"},{"location":"types/#PeriodicGraphEmbeddings.PeriodicGraphEmbedding-Union{Tuple{T}, Tuple{D}, Tuple{PeriodicGraph{D}, AbstractMatrix{T}, Cell}} where {D, T}","page":"Types","title":"PeriodicGraphEmbeddings.PeriodicGraphEmbedding","text":"PeriodicGraphEmbedding{D,T}(graph::PeriodicGraph{D}, placement::AbstractMatrix{T}, cell::Cell=Cell()) where {D,T}\nPeriodicGraphEmbedding{D}(graph::PeriodicGraph{D}, placement::AbstractMatrix{T}, cell::Cell=Cell()) where D\nPeriodicGraphEmbedding(graph::PeriodicGraph{D}, placement::AbstractMatrix{T}, cell::Cell=Cell())\n\nBuild a PeriodicGraphEmbedding{D,T} from the corresponding graph and placement of the vertices, such that each vertex has its fractional coordinate represented in a column of the matrix.\n\nCoordinates out of [0, 1) are translated back to the unit cell with the corresponding offset added to the graph.\n\nThe cell optional argument will not be used if D > 3.\n\nwarning: Warning\nThis function modifies the input graph if any element of placement is out of [0, 1).\n\nnote: Note\nTo obtain a PeriodicGraphEmbedding with sorted positions, use SortedPeriodicGraphEmbedding instead\n\n\n\n\n\n","category":"method"},{"location":"types/#PeriodicGraphEmbeddings.SortedPeriodicGraphEmbedding","page":"Types","title":"PeriodicGraphEmbeddings.SortedPeriodicGraphEmbedding","text":"SortedPeriodicGraphEmbedding{T}\n\nConstructor for PeriodicGraphEmbedding{D,T} where D with sorted positions.\n\n\n\n\n\n","category":"type"},{"location":"types/#PeriodicGraphEmbeddings.SortedPeriodicGraphEmbedding-Union{Tuple{T}, Tuple{D}, Tuple{PeriodicGraph{D}, AbstractMatrix{T} where T, Cell}} where {D, T}","page":"Types","title":"PeriodicGraphEmbeddings.SortedPeriodicGraphEmbedding","text":"SortedPeriodicGraphEmbedding{T}(graph::PeriodicGraph{D}, placement::AbstractMatrix, cell::Cell=Cell()) where {D,T}\n\nBuild a PeriodicGraphEmbedding{D,T} from the corresponding graph and placement of the vertices, so that the result has its vertices sorted by position.\n\nReturn the PeriodicGraphEmbedding as well as the permutation of the columns of placement that yielded the resulting order on the vertices.\n\nThe cell optional argument will not be used if D > 3.\n\nwarning: Warning\nThis function modifies the input graph if any element of placement is out of [0, 1).\n\nSee also PeriodicGraphEmbedding{D,T}(graph, placement::AbstractMatrix{T}, cell) where {D,T} and SortedPeriodicGraphEmbedding(graph, placement::AbstractMatrix, cell).\n\n\n\n\n\nSortedPeriodicGraphEmbedding(graph::PeriodicGraph{D}, placement::AbstractMatrix, cell::Cell=Cell()) where D\n\nBuild a PeriodicGraphEmbedding{D,T} from the corresponding graph and placement of the vertices, so that the result has its vertices sorted by position. T is determined as the smallest type between Rational{Int32}, Rational{Int64}, Rational{Int128} and Rational{BigInt} that can fit all the elements of placement with some additional margin.\n\nReturn the PeriodicGraphEmbedding as well as the permutation of the columns of placement that yielded the resulting order on the vertices.\n\nThe cell optional argument will not be used if D > 3.\n\nwarning: Warning\nThis function modifies the input graph if any element of placement is out of [0, 1).\n\ntip: Tip\nThis function is inherently type-unstable since T cannot be statically determined. This can be useful because having a too large T may slow down later computations.To provide the parameter explicitly, pass it to the SortedPeriodicGraphEmbedding constructor by calling SortedPeriodicGraphEmbedding{T}(graph, placement, cell).\n\nSee also PeriodicGraphEmbedding{D,T}(graph, placement::AbstractMatrix{T}, cell) where {D,T}.\n\n\n\n\n\n","category":"method"},{"location":"types/#PeriodicGraphEmbeddings.PeriodicGraphEmbedding-Union{Tuple{PeriodicGraphEmbedding{N, S}}, Tuple{S}, Tuple{N}, Tuple{T}, Tuple{D}} where {D, T, N, S}","page":"Types","title":"PeriodicGraphEmbeddings.PeriodicGraphEmbedding","text":"PeriodicGraphEmbedding{D,T}(pge::PeriodicGraphEmbedding{N,S}) where {D,T,N,S}\nPeriodicGraphEmbedding{D}(pge::PeriodicGraphEmbedding{N,S}) where {D,N,S}\n\nReturn a PeriodicGraphEmbedding{D,T} with the same structural information as the input pge but embedded in D dimensions instead of N.\n\nIf T is not provided it defaults to S.\n\nThe same caveats that apply to PeriodicGraph{D}(graph::PeriodicGraph{N}) are valid here: namely, the dimensionality of the graph should be at least D and the behaviour is undefined if D < N and there are multiple non-identical connected components.\n\nMoreover, if D < N, the N-D last coordinates of all vertices must be zero or this function will error.\n\n\n\n\n\n","category":"method"},{"location":"types/#Cell-API","page":"Types","title":"Cell API","text":"","category":"section"},{"location":"types/","page":"Types","title":"Types","text":"Cell\ncell_parameters\nEquivalentPosition\nBase.parse(::Type{EquivalentPosition}, s::AbstractString)\nfind_refid","category":"page"},{"location":"types/#PeriodicGraphEmbeddings.Cell","page":"Types","title":"PeriodicGraphEmbeddings.Cell","text":"Cell{T}\n\nRepresentation of a periodic cell in 3D. Contains information about the cell (axes lengths and angles) and its symmetry group, through its Hall number.\n\nSee PeriodicGraphEmbeddings.SPACE_GROUP_HALL, PeriodicGraphEmbeddings.SPACE_GROUP_FULL, PeriodicGraphEmbeddings.SPACE_GROUP_HM and PeriodicGraphEmbeddings.SPACE_GROUP_IT for the correspondance between Hall number and usual symbolic representations.\n\n\n\n\n\n","category":"type"},{"location":"types/#PeriodicGraphEmbeddings.cell_parameters","page":"Types","title":"PeriodicGraphEmbeddings.cell_parameters","text":"cell_parameters(cell::Cell)\n\nReturn ((lengths, angles), mat) where mat is the matrix of the cell in upper triangular format, lengths is the triplet (a, b, c) of lengths of the three axes, and angles is the triplet (α, β, γ) of angles between them.\n\n\n\n\n\n","category":"function"},{"location":"types/#PeriodicGraphEmbeddings.EquivalentPosition","page":"Types","title":"PeriodicGraphEmbeddings.EquivalentPosition","text":"EquivalentPosition{T}\n\nRepresentation of a symmetry operation in 3D, defined by a matrix multiplication and addition.\n\nExample\n\njulia> eq = parse(EquivalentPosition, \"1-x, z, y+1/2\")\n-x+1,z,y+1/2\n\njulia> eq([1//3, 0, 1//4])\n3-element StaticArrays.SVector{3, Rational{Int64}} with indices SOneTo(3):\n 2//3\n 1//4\n 1//2\n\n\n\n\n\n","category":"type"},{"location":"types/#Base.parse-Tuple{Type{EquivalentPosition}, AbstractString}","page":"Types","title":"Base.parse","text":"Base.parse(::Type{EquivalentPosition}, s::AbstractString, refid=(\"x\", \"y\", \"z\"))\n\nParse a string into its represented EquivalentPosition given the name of the three variables obtained from find_refid.\n\nExample\n\njulia> parse(EquivalentPosition, \"a+0.5; c; -1/3+b\", (\"a\", \"b\", \"c\"))\nx+1/2,z,y-1/3\n\njulia> \n\n\n\n\n\n","category":"method"},{"location":"types/#PeriodicGraphEmbeddings.find_refid","page":"Types","title":"PeriodicGraphEmbeddings.find_refid","text":"find_refid(eqs)\n\nFind the reference identifiers for the three dimensions for the CIF group called symmetry_equiv_pos_as_xyz or space_group_symop_operation_xyz. Usually this is simply (\"x\", \"y\", \"z\").\n\n\n\n\n\n","category":"function"},{"location":"utilities/#Utilities","page":"Utilities","title":"Utilities","text":"","category":"section"},{"location":"utilities/#Periodic-distance","page":"Utilities","title":"Periodic distance","text":"","category":"section"},{"location":"utilities/","page":"Utilities","title":"Utilities","text":"The periodic_distance function is useful to compute the shortest distance between two vertices of a graph of dimension 3 or less, i.e. for which a Cell has been provided. It assumes that the unit cell is not too much skewed.","category":"page"},{"location":"utilities/","page":"Utilities","title":"Utilities","text":"periodic_distance","category":"page"},{"location":"utilities/#PeriodicGraphEmbeddings.periodic_distance","page":"Utilities","title":"PeriodicGraphEmbeddings.periodic_distance","text":"periodic_distance(u, mat, ortho=nothing, safemin=nothing)\n\nDistance between point u and the origin, given as a triplet of fractional coordinates, in a repeating unit cell of matrix mat. The distance is the shortest between all equivalents of u and the origin. If ortho is set to true, the angles α, β and γ of the cell are assumed right, which accelerates the computation by up to 7 times. If a distance lower than safemin is computed, stop trying to find a periodic image of u closer to the origin. If unspecified, both ortho and safemin are automatically determined from mat.\n\nThis implementation assumes that the cell corresponds to a reduced lattice. It may be invalid for some edge cases otherwise.\n\nFor optimal performance, use periodic_distance! with buffer, ortho and safemin obtained from prepare_periodic_distance_computations.\n\n\n\n\n\n","category":"function"},{"location":"utilities/#Other","page":"Utilities","title":"Other","text":"","category":"section"},{"location":"utilities/","page":"Utilities","title":"Utilities","text":"PeriodicGraphEmbeddings.double_widen","category":"page"},{"location":"utilities/#PeriodicGraphEmbeddings.double_widen","page":"Utilities","title":"PeriodicGraphEmbeddings.double_widen","text":"double_widen(::Type)\n\nInternal function used to selectively widen small integer and rational types.\n\nThis is useful to avoid overflow without sacrificing too much efficiency by always having to resolve to very large types.\n\n\n\n\n\n","category":"function"},{"location":"#PeriodicGraphEmbeddings","page":"Home","title":"PeriodicGraphEmbeddings","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A Julia package for manipulating periodic graph embeddings in euclidean spaces, such as those representing atoms in a crystal.","category":"page"},{"location":"","page":"Home","title":"Home","text":"See also:","category":"page"},{"location":"","page":"Home","title":"Home","text":"PeriodicGraphs.jl for the all that relates to the periodic graph itself, irrespective of its euclidean embedding.\nCrystalNets.jl for a dependent package specialized on crystal nets.","category":"page"},{"location":"io/#I/O","page":"I/O","title":"I/O","text":"","category":"section"},{"location":"io/","page":"I/O","title":"I/O","text":"export_vtf\nexport_cgd","category":"page"},{"location":"io/#PeriodicGraphEmbeddings.export_vtf","page":"I/O","title":"PeriodicGraphEmbeddings.export_vtf","text":"function export_vtf(file, pge::PeriodicGraphEmbedding3D{T}, types=nothing, repeatedges=6, colorname=false, tostring=string, atomnumof==(a,i)->(a isa Integer ? a : i)) where T\n\nExport a PeriodicGraphEmbedding3D to a .vtf file (readable by VMD).\n\nIf specified, types is a list of types for each vertex of pge. Each type is converted to string by the tostring function. The atomnumof function takes two arguments ty and i where ty is a type and i is the number of the vertex, and return an Int representing an atom number.\n\n\n\n\n\n","category":"function"},{"location":"io/#PeriodicGraphEmbeddings.export_cgd","page":"I/O","title":"PeriodicGraphEmbeddings.export_cgd","text":"export_cgd(file, pge::PeriodicGraphEmbedding, name=basename(splitext(file)[1]), append=false)\nexport_cgd(file, g::PeriodicGraph, name=basename(splitext(file)[1]), append=false)\n\nExport a PeriodicGraph or a PeriodicGraphEmbedding to a .cgd file (readable by Systre).\n\nIf append is set, the graph is added at the end of the file.\n\n\n\n\n\n","category":"function"},{"location":"symmetries/#Symmetries","page":"Symmetries","title":"Symmetries","text":"","category":"section"},{"location":"symmetries/","page":"Symmetries","title":"Symmetries","text":"Symmetry detection is provided by the spglib library, wrapped in helper functions detailed below.","category":"page"},{"location":"symmetries/#Manual","page":"Symmetries","title":"Manual","text":"","category":"section"},{"location":"symmetries/","page":"Symmetries","title":"Symmetries","text":"The main function is find_symmetries which returns a SymmetryGroup3D:","category":"page"},{"location":"symmetries/","page":"Symmetries","title":"Symmetries","text":"find_symmetries\nSymmetryGroup3D\nPeriodicGraphEmbeddings.PeriodicSymmetry3D","category":"page"},{"location":"symmetries/#PeriodicGraphEmbeddings.find_symmetries","page":"Symmetries","title":"PeriodicGraphEmbeddings.find_symmetries","text":"find_symmetries(pge::PeriodicGraphEmbedding3D, vtypes=nothing, check_symmetry=check_valid_symmetry)\n\nReturn a SymmetryGroup3D object storing the list of symmetry operations on the graph embedding.\n\nIf vtypes !== nothing, ensure that two vertices x and y cannot be symmetry-related if vtypes[x] != vtypes[y].\n\ncheck_symmetry must be a function that takes the same four arguments pge, t, r and vtypes as check_valid_symmetry and return either (vmap, offsets) or nothing if the input is not a valid symmetry. It can be used to specify additional constraints that cannot be carried by vtypes alone.\n\n\n\n\n\n","category":"function"},{"location":"symmetries/#PeriodicGraphEmbeddings.SymmetryGroup3D","page":"Symmetries","title":"PeriodicGraphEmbeddings.SymmetryGroup3D","text":"SymmetryGroup3D{T} <: PeriodicGraphs.AbstractSymmetryGroup\n\nStore the information on the symmetry operations available on a PeriodicGraphEmbedding3D.\n\n\n\n\n\n","category":"type"},{"location":"symmetries/#PeriodicGraphEmbeddings.PeriodicSymmetry3D","page":"Symmetries","title":"PeriodicGraphEmbeddings.PeriodicSymmetry3D","text":"PeriodicSymmetry3D{T} <: PeriodicGraphs.AbstractSymmetry\n\nSingle symmetry of a PeriodicGraphEmbedding3D{T}.\n\nSee PeriodicGraphs.AbstractSymmetry for information on the API.\n\n\n\n\n\n","category":"type"},{"location":"symmetries/#Space-group-database-API","page":"Symmetries","title":"Space group database API","text":"","category":"section"},{"location":"symmetries/","page":"Symmetries","title":"Symmetries","text":"PeriodicGraphEmbeddings.SPACE_GROUP_HALL\nPeriodicGraphEmbeddings.SPACE_GROUP_HM\nPeriodicGraphEmbeddings.SPACE_GROUP_FULL\nPeriodicGraphEmbeddings.SPACE_GROUP_IT\nPeriodicGraphEmbeddings.HALL_SYMBOLS","category":"page"},{"location":"symmetries/#PeriodicGraphEmbeddings.SPACE_GROUP_HALL","page":"Symmetries","title":"PeriodicGraphEmbeddings.SPACE_GROUP_HALL","text":"Dictionnary mapping the Hall symbol of a symmetry group to its Hall number.\n\nIn the keys, letters are lowercase, underscores are removed and space is kept to differentiate \"p 6 2\" from \"p 62\" and \"p 3 2\" from \"p 32\"\n\n\n\n\n\n","category":"constant"},{"location":"symmetries/#PeriodicGraphEmbeddings.SPACE_GROUP_HM","page":"Symmetries","title":"PeriodicGraphEmbeddings.SPACE_GROUP_HM","text":"Dictionnary mapping the HM symbol of a symmetry group to its Hall number.\n\nIn the keys, letters are lowercase and space is removed.\n\n\n\n\n\n","category":"constant"},{"location":"symmetries/#PeriodicGraphEmbeddings.SPACE_GROUP_FULL","page":"Symmetries","title":"PeriodicGraphEmbeddings.SPACE_GROUP_FULL","text":"Dictionnary mapping the full notation representation of a symmetry group to its Hall number, if the full notation is distinct from the H-M symbol.\n\nIn the keys, letters are lowercase and space is removed.\n\n\n\n\n\n","category":"constant"},{"location":"symmetries/#PeriodicGraphEmbeddings.SPACE_GROUP_IT","page":"Symmetries","title":"PeriodicGraphEmbeddings.SPACE_GROUP_IT","text":"List mapping the International Table number of a symmetry group to its Hall number\n\n\n\n\n\n","category":"constant"},{"location":"symmetries/#PeriodicGraphEmbeddings.HALL_SYMBOLS","page":"Symmetries","title":"PeriodicGraphEmbeddings.HALL_SYMBOLS","text":"List of Hall symbols and crystal system corresponding to each Hall number\n\n\n\n\n\n","category":"constant"},{"location":"symmetries/#Internal-API","page":"Symmetries","title":"Internal API","text":"","category":"section"},{"location":"symmetries/","page":"Symmetries","title":"Symmetries","text":"find_hall_number\nPeriodicGraphEmbeddings.SpglibDataset\nget_symmetry_equivalents\nPeriodicGraphEmbeddings.get_spglib_dataset\ncheck_valid_symmetry","category":"page"},{"location":"symmetries/#PeriodicGraphEmbeddings.find_hall_number","page":"Symmetries","title":"PeriodicGraphEmbeddings.find_hall_number","text":"find_hall_number(hallsymbol::AbstractString, hm::AbstractString=hallsymbol, it::Integer=0, warnonnotfound=false)\n\nDetermine the hall number corresponding to the given hallsymbol. The Hermann-Mauguin symbol hm can alternatively be used, or simply the International Table number of the space group it to get the hall number of the standard setting of the group.\n\nPassing an empty string to hallsymbol or hm or 0 to it disregards the argument.\n\nThe optional argument warnonnotfound specifies whether to print a warning if one of the provided arguments was not reckognized.\n\n\n\n\n\n","category":"function"},{"location":"symmetries/#PeriodicGraphEmbeddings.SpglibDataset","page":"Symmetries","title":"PeriodicGraphEmbeddings.SpglibDataset","text":"SpglibDataset\n\nWrapper around the SpglibDataset type exported by spglib. Its accessible fields are the same as in the C counterpart, except that strings are already converted to String, lists to Vector and matrices to Matrix.\n\nTo access the raw pointers without conversion, prepend an underscore to the field: for example dataset._rotations yields a Ptr{Cint} where dataset.rotations is a 3×3 Matrix{Int}.\n\n\n\n\n\n","category":"type"},{"location":"symmetries/#PeriodicGraphEmbeddings.get_symmetry_equivalents","page":"Symmetries","title":"PeriodicGraphEmbeddings.get_symmetry_equivalents","text":"get_symmetry_equivalents([T=Rational{Int},] hall)\n\nThe list of EquivalentPosition{T} corresponding to a symmetry group given by its Hall number.\n\nWrapper around spg_get_symmetry_from_database.\n\n\n\n\n\n","category":"function"},{"location":"symmetries/#PeriodicGraphEmbeddings.get_spglib_dataset","page":"Symmetries","title":"PeriodicGraphEmbeddings.get_spglib_dataset","text":"get_spglib_dataset(pge::PeriodicGraphEmbedding3D, vtypes=nothing)\n\nWrapper around spg_get_dataset.\n\nIf vtypes !== nothing, ensure that two vertices x and y cannot be symmetry-related if vtypes[x] != vtypes[y].\n\n\n\n\n\n","category":"function"},{"location":"symmetries/#PeriodicGraphEmbeddings.check_valid_symmetry","page":"Symmetries","title":"PeriodicGraphEmbeddings.check_valid_symmetry","text":"check_valid_symmetry(pge::PeriodicGraphEmbedding{D,T}, t::SVector{D,T}, r=nothing, vtypes=nothing, issorted=false)\n\nCheck that the periodic graph embedding is identical to that rotated by r (if it is not nothing) then translated by t. If vtypes is not nothing, any vertex x must additionally be mapped to a vertex y such that vtypes[x] == vtypes[y]. If issorted is set and T <: Rational, assume that issorted(pge.pos) to use a faster dichotomy approach.\n\nIf so, return the the vmap between the initial vertices and their symmetric images, as well as the offsets of each symmetric image compared to the origin. Otherwise, return nothing.\n\n\n\n\n\n","category":"function"}]
}