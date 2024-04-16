# Other utilies from PeriodicGraphs.jl adapted to PeriodicGraphEmbeddings

function PeriodicGraphs.offset_representatives!(pge::PeriodicGraphEmbedding, offsets)
    offset_representatives!(pge.g, offsets)
    for (i,x) in enumerate(pge.pos)
        pge.pos[i] = x .- offsets[i]
    end
    nothing
end

function PeriodicGraphs.swap_axes!(pge::PeriodicGraphEmbedding{N}, t) where N
    swap_axes!(pge.g, t)
    for (i,x) in enumerate(pge.pos)
        pge.pos[i] = x[t]
    end
    nothing
end

function PeriodicGraphs.make_supercell(cell::Cell{T}, t) where T
    newmat, newequivalents = if length(t) == 2
        isempty(cell.equivalents) || error("2D supercell including symmetries are not currently handled")
        SMatrix{3,3,BigFloat,9}(
            cell.mat[1,1]*t[1], cell.mat[2,1]*t[1], cell.mat[3,1]*t[1],
            cell.mat[1,2]*t[2], cell.mat[2,2]*t[2], cell.mat[3,2]*t[2],
            cell.mat[1,3], cell.mat[2,3], cell.mat[3,3],
        ), copy(cell.equivalents)
    elseif length(t) == 3
        t1, t2, t3 = t
        newmat = SMatrix{3,3,BigFloat,9}(
            cell.mat[1,1]*t1, cell.mat[2,1]*t1, cell.mat[3,1]*t1,
            cell.mat[1,2]*t2, cell.mat[2,2]*t2, cell.mat[3,2]*t2,
            cell.mat[1,3]*t3, cell.mat[2,3]*t3, cell.mat[3,3]*t3,
        )
        n = length(cell.equivalents)
        if n == 0
            newmat, EquivalentPosition{T}[]
        else
            newequivalents = Vector{EquivalentPosition{T}}(undef, n + (n+1)*(prod(t)-1))
            it1, it2, it3 = 1//t1, 1//t2, 1//t3
            for (i, eq) in enumerate(cell.equivalents)
                newequivalents[i] = EquivalentPosition(SMatrix{3,3,T,9}(
                    eq.mat[1,1],        eq.mat[2,1]*t1*it2, eq.mat[3,1]*t1*it3,
                    eq.mat[1,2]*t2*it1, eq.mat[2,2],        eq.mat[3,2]*t2*it3,
                    eq.mat[1,3]*t3*it1, eq.mat[2,3]*t3*it2, eq.mat[3,3]
                ), SVector{3,T}(eq.ofs[1]*it1, eq.ofs[2]*it2, eq.ofs[3]*it3))
            end
            for (j, (a, b, c)) in enumerate(PeriodicGraphs.MetaClock(t))
                a == b == c == 0 && continue
                ofs = (j-1)*(n+1)
                trans = SVector{3,T}(a*it1, b*it2, c*it3)
                newequivalents[ofs] = EquivalentPosition(one(SMatrix{3,3,T,9}), trans)
                for i in 1:n
                    ref = newequivalents[i]
                    newequivalents[ofs+i] = EquivalentPosition(ref.mat, ref.ofs + trans)
                end
            end
            newmat, newequivalents
        end
    end
    Cell{T}(cell.hall, newmat, newequivalents)
end
function PeriodicGraphs.make_supercell(pge::PeriodicGraphEmbedding, t)
    g = make_supercell(pge.g, t)
    cell = make_supercell(pge.cell, t)
    newpos = copy(pge.pos)
    n = length(newpos)
    for i in 1:n
        newpos[i] = newpos[i] ./ t
    end
    clock = PeriodicGraphs.MetaClock(t)
    resize!(newpos, n*length(clock))
    for (k, pos) in enumerate(Iterators.drop(clock, 1))
        for i in 1:n
            newpos[n*k+i] = newpos[i] .+ pos .// t
        end
    end
    return PeriodicGraphEmbedding(g, newpos, cell)
end
