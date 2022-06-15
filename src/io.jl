export export_vtf, export_cgd

"""
    function export_vtf(file::AbstractString, pge::PeriodicGraphEmbedding3D{T}, types=nothing, repeatedges=6, colorname=false, tostring=string, atomnumof==(a,i)->(a isa Integer ? a : i)) where T

Export a [`PeriodicGraphEmbedding3D`](@ref) to a .vtf `file` (readable by VMD).

If specified, `types` is a list of types for each vertex of `pge`. Each type is converted
to string by the `tostring` function.
The `atomnumof` function takes two arguments `ty` and `i` where `ty` is a type and `i` is
the number of the vertex, and return an `Int` representing an atom number.
"""
function export_vtf(file::AbstractString, pge::PeriodicGraphEmbedding3D{T}, types=nothing, repeatedges=6, colorname=false, tostring=string, atomnumof=(a,i)->(a isa Integer ? a : i)) where T
    endswith(file, ".vtf") || return export_vtf(file*".vtf", pge, types, repeatedges, colorname, tostring, atomnumof)
    mkpath(splitdir(file)[1])
    n = nv(pge.g)
    open(file, write=true) do f
        invcorres = [PeriodicVertex3D(i) for i in 1:n]
        corres = Dict{PeriodicVertex3D,Int}([invcorres[i]=>i for i in 1:n])
        encounteredtypes = Dict{Symbol,String}()
        atomnums = Int[]
        numencounteredtypes = 0
        println(f, """
        ###############################
        # written by PeriodicGraphEmbeddings.jl
        ###############################
        """)

        for i in 1:n
            ty = isnothing(types) ? Symbol("") : types[i]
            sty = ty === Symbol("") ? string(i) : tostring(ty)
            if length(sty) > 16
                sty = sty[1:13]*"etc" # otherwise VMD fails to load the .vtf
            end
            name = if colorname
                _name = get(encounteredtypes, ty, missing)
                if _name isa String
                    _name
                else
                    numencounteredtypes += 1
                    encounteredtypes[ty] = string(" name ", numencounteredtypes)
                end
            else
                n ≥ 32768 ? "" : string(" name ", i)
            end
            atomnum = ty === Symbol("") ? 0 : atomnumof(ty, i)
            push!(atomnums, atomnum)
            resid = colorname ? i : 0
            println(f, "atom $(i-1) type $sty$name resid $resid atomicnumber $atomnum")
        end
        j = n + 1
        for _ in 1:repeatedges
            jmax = j - 1
            for i in 1:jmax
                vertex = invcorres[i]
                for x in neighbors(pge.g, vertex.v)
                    y = PeriodicVertex3D(x.v, x.ofs .+ vertex.ofs)
                    if get!(corres, y, j) == j
                        j += 1
                        push!(invcorres, y)
                    end
                end
            end
        end
        for i in n+1:length(invcorres)
            v = invcorres[i].v
            ofs = invcorres[i].ofs
            ty = isnothing(types) ? Symbol("") : types[v]
            sty = ty === Symbol("") ? string(i) : tostring(ty)
            if length(sty) > 16
                sty = sty[1:13]*"etc"
            end
            name = colorname ? string(" name ", encounteredtypes[ty]) :
                    n ≥ 32768 ? "" : string(" name ", v)
            atomnum = atomnums[v]
            resid = colorname ? i : PeriodicGraphs.hash_position(ofs)
            println(f, "atom $(i-1) type $sty$name resid $resid atomicnumber $atomnum")
        end
        println(f)

        for (i,x) in enumerate(invcorres)
            for neigh in neighbors(pge.g, x.v)
                j = get(corres, PeriodicVertex3D(neigh.v, neigh.ofs .+ x.ofs), nothing)
                isnothing(j) && continue
                if i < j
                    println(f, "bond ", i - 1, ':', j - 1)
                end
            end
        end
        println(f)

        ((_a, _b, _c), (_α, _β, _γ)), mat = cell_parameters(pge.cell)
        println(f, "pbc $_a $_b $_c $_α $_β $_γ\n")

        println(f, "ordered")
        for x in invcorres
            coord = mat * ((T <: Rational ? widen.(pge.pos[x.v]) : pge.pos[x.v]) .+ x.ofs)
            join(f, round.(Float64.(coord); digits=15), ' ')
            println(f)
        end
    end
    nothing
end

"""
    export_cgd(file, pge::PeriodicGraphEmbedding, name=basename(splitext(file)[1]), append=false)
    export_cgd(file, g::PeriodicGraph, name=basename(splitext(file)[1]), append=false)

Export a `PeriodicGraph` or a `PeriodicGraphEmbedding` to a .cgd `file` (readable by
Systre).

If `append` is set, the graph is added at the end of the file.
"""
function export_cgd(file, pge::PeriodicGraphEmbedding, name=basename(splitext(file)[1]), append=false)
    endswith(file, ".cgd") || return export_cgd(file*".cgd", pge, name, append)
    mkpath(splitdir(file)[1])
    open(file; write=true, append) do f
        println(f, "CRYSTAL")
        println(f, "  NAME\t", name)
        ((__a, __b, __c), (__α, __β, __γ)), _ = cell_parameters(pge.cell)
        _a, _b, _c, _α, _β, _γ = Float64.((__a, __b, __c, __α, __β, __γ))
        print(f, "  GROUP\t\"")
        print(f, RAW_SYMMETRY_DATA[pge.cell.hall][4])
        println(f, "\"")
        println(f, "  CELL\t", _a, ' ', _b, ' ', _c, ' ', _α, ' ', _β, ' ', _γ, ' ')
        println(f, "  ATOM")
        for i in 1:length(pge.pos)
            pos = pge.pos[i]
            println(f, "    ", i, ' ', degree(pge.g, i), ' ', pos[1], ' ',
                    pos[2], ' ', pos[3])
        end
        println(f, "  EDGE")
        for i in 1:length(pge.pos)
            for e in neighbors(pge.g, i)
                e.v < i && continue
                dest = pge.pos[e.v] .+ e.ofs
                println(f, "    ", i, '\t', dest[1], ' ', dest[2], ' ', dest[3])
            end
        end
        println(f, "\nEND")
    end
end

function export_cgd(file, g::PeriodicGraph, name=basename(splitext(file)[1]), append=false)
    endswith(file, ".cgd") || return export_cgd(file*".cgd", pge, name, append)
    mkpath(splitdir(file)[1])
    open(file; write=true, append) do f
        println(f, "PERIODIC_GRAPH")
        println(f, "  NAME ", name)
        println(f, "  EDGES")
        repr = reverse(split(string(g)))
        n = parse(Int, pop!(repr))
        m = length(repr) ÷ (n+2)
        for _ in 1:m
            src = pop!(repr)
            dst = pop!(repr)
            ofs = Vector{String}(undef, n)
            for i in 1:n
                ofs[i] = pop!(repr)
            end
            print(f, "    ", src, ' ', dst, ' ')
            join(f, ofs, ' ')
            println(f)
        end
        println(f, "END\n")
    end
end
