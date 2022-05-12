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
end