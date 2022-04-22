export prepare_periodic_distance_computations, periodic_distance!, periodic_distance

@static if VERSION < v"1.8-"
    # copy-pasted from the implementation of @lazy_str in Base
    macro lazy_str(text)
        parts = Any[]
        lastidx = idx = 1
        while (idx = findnext('$', text, idx)) !== nothing
            lastidx < idx && push!(parts, text[lastidx:idx-1])
            idx += 1
            expr, idx = Meta.parseatom(text, idx; filename=string(__source__.file))
            push!(parts, esc(expr))
            lastidx = idx
        end
        lastidx <= lastindex(text) && push!(parts, text[lastidx:end])
        :(string($(parts...)))
    end
end



function prepare_periodic_distance_computations(mat)
    (a, b, c), (α, β, γ) = cell_parameters(mat)
    ortho = all(x -> isapprox(Float16(x), 90; rtol=0.02), (α, β, γ))
    _a, _b, _c = eachcol(mat)
    safemin = min(Float64(dot(cross(_b, _c), _a)/(b*c)),
                  Float64(dot(cross(_c, _a), _b)/(a*c)),
                  Float64(dot(cross(_a, _b), _c)/(a*b)))/2
    # safemin is the half-distance between opposite planes of the unit cell
    return MVector{3,Float64}(undef), ortho, safemin
end

function periodic_distance!(buffer, u, mat, ortho, safemin)
    @simd for i in 1:3
        diff = u[i] + 0.5
        buffer[i] = diff - floor(diff) - 0.5
    end
    ref = norm(mat*buffer)
    (ortho || ref ≤ safemin) && return ref
    @inbounds for i in 1:3
        buffer[i] += 1
        newnorm = norm(mat*buffer)
        newnorm < ref && return newnorm # in a reduced lattice, there should be at most one
        buffer[i] -= 2
        newnorm = norm(mat*buffer)
        newnorm < ref && return newnorm
        buffer[i] += 1
    end
    return ref
end

"""
    periodic_distance(u, mat, ortho=nothing, safemin=nothing)

Distance between point `u` and the origin, given as a triplet of fractional coordinates, in
a repeating unit cell of matrix `mat`.
The distance is the shortest between all equivalents of `u` and the origin.
If `ortho` is set to `true`, the angles α, β and γ of the cell are assumed right, which
accelerates the computation by up to 7 times.
If a distance lower than `safemin` is computed, stop trying to find a periodic image of `u`
closer to the origin.
If unspecified, both `ortho` and `safemin` are automatically determined from `mat`.

This implementation assumes that the cell corresponds to a reduced lattice. It may be
invalid for some edge cases otherwise.

For optimal performance, use `periodic_distance!` with `buffer`, `ortho` and `safemin`
obtained from `prepare_periodic_distance_computations`.
"""
function periodic_distance(u, mat, ortho=nothing, safemin=nothing)
    if ortho === nothing || safemin === nothing
        _, ortho, safemin = prepare_periodic_distance_computations(mat)
    end
    periodic_distance!(similar(u), u, mat, ortho::Bool, safemin::Float64)
end

