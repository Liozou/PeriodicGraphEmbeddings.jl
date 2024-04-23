export PeriodicDistance2
using LinearAlgebra: mul!

@static if VERSION < v"1.7-"
    Base.replace(x::String, p1::Pair, p2::Pair...) = replace(replace(x, p1), p2...)
end

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

function norm2(u::T) where {T}
    r2 = zero(eltype(T))^2
    @simd for x in u
        r2 += x^2
    end
    r2
end

"""
    (pd2::PeriodicDistance2)(x, y=nothing, ofsx=nothing, ofsy=nothing; fromcartesian=false, ofs=nothing)

Squared periodic distance between points `x` and `y`, optionally with respective offsets
`ofsx` and `ofsy`.
The periodic distance is the shortest distance between all the periodic images of the
inputs. For `x` and `y` given as triplets of fractional coordinates, putting integer
offsets in `ofsx` and `ofsy` does not have any effect thus.

If `y` is not set, compute the periodic distance between `x` and the origin.

If `fromcartesian` is set, the inputs must be given as cartesian coordinates. Otherwise, it
is assumed that the input is given in fractional coordinates.

If `ofs` is set to a mutable vector of integers (`MVector{3,Int}` from StaticArrays.jl is
suggested), then the offset of the translation of `y .+ ofsy` with respect to `x .+ ofsx`
is stored in `ofs`.
This means that the returned periodic distance is equal to the non-periodic distance
between `y .+ ofsy .+ ofs` and `x .+ ofsx` (assuming fractional coordinates).
The offset is always returned as a triplet of integers, though the input can be given in
cartesian coordinates as long as `fromcartesian` is set.

!!! warning
    This function modifies an internal  state of `pd2`, and is thus not thread-safe.

## Example
```jldoctest
julia> mat = [26.04 -7.71 -8.32; 0.0 32.72 -3.38; 0.0 0.0 26.66];

julia> pd2 = PeriodicDistance2(mat)
PeriodicDistance2([26.04 -7.71 -8.32; 0.0 32.72 -3.38; 0.0 0.0 26.66])

julia> sqrt(pd2([0.5, 0, 0])) # distance to the middle of the a axis is simply a/2
13.02

julia> ofs = MVector{3,Int}(undef);

julia> vec1 = [0.9, 0.6, 0.5]; vec2 = [0.0, 0.5, 0.01];

julia> d2 = pd2(vec1, vec2; ofs)
210.57932044000003

julia> println(ofs) # the offset of vec2 to have the closest image to vec1
[1, 0, 1]

julia> norm(mat*(vec2 .+ ofs .- vec1))^2 ≈ d2
true

julia> pd2(mat*vec1, mat*vec2; fromcartesian=true) ≈ d2
true
```
"""
struct PeriodicDistance2{T,Ti,T2}
    buffer::MVector{3,T}
    buffer2::MVector{3,Float64}
    mat::SMatrix{3,3,T,9}
    invmat::SMatrix{3,3,Ti,9}
    ortho::Bool
    safemin2::T2
end
Base.show(io::IO, pd2::PeriodicDistance2) = println(io, PeriodicDistance2, '(', pd2.mat, ')')

"""
    PeriodicDistance2(mat::AbstractMatrix)

Build a [`PeriodicDistance2`](@ref) object which can be called to compute squared periodic
distances in a unit cell of given matrix `mat`.

!!! warning
    It is assumed that the unit cell is in standard settings: its angles must be above 60°
    or the returned distance may not be the shortest.
"""
function PeriodicDistance2(mat::AbstractMatrix{T}) where {T}
    smat = SMatrix{3,3,T,9}(mat)
    (a, b, c), (α, β, γ) = cell_parameters(smat)
    ortho = all(x -> isapprox(Float16(x), 90; rtol=0.02), (α, β, γ))
    _a, _b, _c = eachcol(smat)
    safemin = min(Float64(dot(cross(_b, _c), _a)/(b*c)),
                  Float64(dot(cross(_c, _a), _b)/(a*c)),
                  Float64(dot(cross(_a, _b), _c)/(a*b)))/2
    # safemin is the half-distance between opposite planes of the unit cell
    buffer = MVector{3,T}(undef)
    buffer2 = MVector{3,Float64}(undef)
    T1 = oneunit(first(smat))
    invmat = inv(smat./T1)./T1 # to handle Unitful
    PeriodicDistance2(buffer, buffer2, smat, invmat, ortho, safemin^2)
end

@inline function _set_if_ofs!(buffer, x, y, ofsx, ofsy)
    if y isa Nothing
        if ofsx isa Nothing
            if ofsy isa Nothing
                buffer .= x
            else
                buffer .= x .- ofsy
            end
        else
            if ofsy isa Nothing
                buffer .= x .+ ofsx
            else
                buffer .= x .+ ofsx .- ofsy
            end
        end
    else
        if ofsx isa Nothing
            if ofsy isa Nothing
                buffer .= x .- y
            else
                buffer .= x .- y .- ofsy
            end
        else
            if ofsy isa Nothing
                buffer .= x .+ ofsx .- y
            else
                buffer .= x .+ ofsx .- y .- ofsy
            end
        end
    end
end

function _periodic_distance2!(pd2::PeriodicDistance2, ofs)
    if ofs isa Nothing
        @simd for i in 1:3
            diff = pd2.buffer2[i] + 0.5
            pd2.buffer2[i] = diff - floor(diff) - 0.5
        end
    else
        @simd for i in 1:3
            diff = pd2.buffer2[i] + 0.5
            tmp = floor(Int, diff)
            ofs[i] = tmp
            pd2.buffer2[i] = diff - tmp - 0.5
        end
    end
    mul!(pd2.buffer, pd2.mat, pd2.buffer2)
    ref2 = norm2(pd2.buffer)
    (pd2.ortho || ref2 ≤ pd2.safemin2) && return ref2
    @inbounds for i in 1:3
        pd2.buffer2[i] += 1
        mul!(pd2.buffer, pd2.mat, pd2.buffer2)
        newnorm2 = norm2(pd2.buffer)
        if newnorm2 < ref2
            ofs isa Nothing || (ofs[i] -= 1)
            return newnorm2 # in a reduced lattice, there should be at most one
        end
        pd2.buffer2[i] -= 2
        mul!(pd2.buffer, pd2.mat, pd2.buffer2)
        newnorm2 = norm2(pd2.buffer)
        if newnorm2 < ref2
            ofs isa Nothing || (ofs[i] += 1)
            return newnorm2
        end
        pd2.buffer2[i] += 1
    end
    return ref2
end

function (pd2::PeriodicDistance2)(x, y=nothing, ofsx=nothing, ofsy=nothing; fromcartesian=false, ofs=nothing)
    if fromcartesian
        _set_if_ofs!(pd2.buffer, x, y, ofsx, ofsy)
        mul!(pd2.buffer2, pd2.invmat, pd2.buffer)
    else
        _set_if_ofs!(pd2.buffer2, x, y, ofsx, ofsy)
    end
    _periodic_distance2!(pd2, ofs)
end
