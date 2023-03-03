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
    newmat = SizedMatrix{3,3,BigFloat,2,Matrix{BigFloat}}(cell.mat)
    for (i,ti) in zip(1:3, t)
        newmat[:,i] .*= ti
    end
    Cell{T}(cell.hall, newmat, cell.equivalents)
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
