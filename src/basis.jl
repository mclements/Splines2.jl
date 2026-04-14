# the original Splines2 implementation used zero-indexed OffsetArrays
# to ease the port from C
struct SplineBasis{T<:Real} <: AbstractSplineBasis{T}
    order::Int # order of the spline
    nknots::Int # number of knots
    ncoef::Int # number of coefficients
    knots::OffsetArray{T,1} # knot vector
end

function SplineBasis(knots::AbstractVector{T}, order::Int=4) where {T<:Real}
    return SplineBasis(order,
                       length(knots),
                       length(knots) - order,
                       zeroIndexedArray(knots))
end

function basis(bs::SplineBasis{T}, x::T, derivs::Int=0) where {T<:Real}
    t = bs.knots
    k = bs.order - 1
    m = derivs
    hh = bs.order
    ell = _find_interval(bs, x)
    result = zeroIndexedArray(zeros(T, 2 * k + 2))
    one = T(1)
    zero = T(0)
    result[0] = one
    for j in 1:(k - m)
        for n in 0:(j - 1)
            result[hh + n] = result[n]
        end
        result[0] = zero
        for n in 1:j
            ind = ell + n
            xb = t[ind]
            xa = t[ind - j]
            if (xb == xa)
                result[n] = zero
                continue
            end
            w = result[hh + n - 1] / (xb - xa)
            result[n - 1] = result[n - 1] + w * (xb - x)
            result[n] = w * (x - xa)
        end
    end
    for j in (k - m + 1):k
        for n in 0:(j - 1)
            result[hh + n] = result[n]
        end
        result[0] = zero
        for n in 1:j
            ind = ell + n
            xb = t[ind]
            xa = t[ind - j]
            if xb == xa
                result[m] = zero
                continue
            end
            w = j * result[hh + n - 1] / (xb - xa)
            result[n - 1] = result[n - 1] - w
            result[n] = w
        end
    end
    offset = ell - k
    v = zeros(T, bs.ncoef)
    v[(1 + offset):(k + 1 + offset)] = result[0:k]
    return v
end

function basis(bs::AbstractSplineBasis{T}, x::AbstractVector{T},
               derivs::Int=0) where {T<:Real}
    f(xi) = basis(bs, xi, derivs)
    return stack(f, x; dims=1)
end
