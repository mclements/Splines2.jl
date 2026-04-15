struct SplineBasis{T<:Real} <: AbstractSplineBasis{T}
    order::Int # order of the spline
    nknots::Int # number of knots
    ncoef::Int # number of coefficients
    knots::Vector{T} # knot vector
end

function SplineBasis(knots::AbstractVector{T}, order::Int=4) where {T<:Real}
    return SplineBasis(order,
                       length(knots),
                       length(knots) - order,
                       Vector{T}(knots))
end

function basis(bs::SplineBasis{T}, x::T, derivs::Int=0) where {T<:Real}
    t = bs.knots
    k = bs.order - 1
    m = derivs
    hh = bs.order
    ell = _find_interval(bs, x)
    result = zeros(T, 2 * k + 2)
    one = T(1)
    zero = T(0)
    result[1] = one
    for j in 1:(k - m)
        for n in 0:(j - 1)
            result[hh + n + 1] = result[n + 1]
        end
        result[1] = zero
        for n in 1:j
            ind = ell + n
            xb = t[ind]
            xa = t[ind - j]
            if (xb == xa)
                result[n + 1] = zero
                continue
            end
            w = result[hh + n] / (xb - xa)
            result[n] = result[n] + w * (xb - x)
            result[n + 1] = w * (x - xa)
        end
    end
    for j in (k - m + 1):k
        for n in 0:(j - 1)
            result[hh + n + 1] = result[n + 1]
        end
        result[1] = zero
        for n in 1:j
            ind = ell + n
            xb = t[ind]
            xa = t[ind - j]
            if xb == xa
                result[m + 1] = zero
                continue
            end
            w = j * result[hh + n] / (xb - xa)
            result[n] = result[n] - w
            result[n + 1] = w
        end
    end
    v = zeros(T, bs.ncoef)
    copyto!(view(v, (ell - k):(ell)), view(result, 1:(k + 1)))
    return v
end

function basis(bs::AbstractSplineBasis{T}, x::AbstractVector{T},
               derivs::Int=0) where {T<:Real}
    f(xi) = basis(bs, xi, derivs)
    return stack(f, x; dims=1)
end
