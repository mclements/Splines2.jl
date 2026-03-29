struct MSplineBasis{T<:Real} <: AbstractSplineBasis{T}
    b_spline_basis::BSplineBasis{T}
    trans_coef::Vector{T}
end

function MSplineBasis(boundary_knots::Tuple{T,T},
                      interior_knots::Union{Array{T,1},Nothing}=nothing,
                      order::Int=4,
                      intercept::Bool=false) where {T<:Real}
    spline = BSplineBasis(boundary_knots, interior_knots, order, intercept)
    knots = spline.spline_basis.knots
    trans_coef = map(1:(length(knots) - order)) do j
        denom = knots[j + order] - knots[j]
        return denom > T(0) ? order / denom : T(0)
    end
    if !intercept
        trans_coef = trans_coef[2:length(trans_coef)]
    end
    return MSplineBasis(spline, trans_coef)
end

function basis(ms::MSplineBasis{T}, x::T, derivs::Int=0) where {T<:Real}
    vec = basis(ms.b_spline_basis, x, derivs)
    vec .*= ms.trans_coef
    return vec
end

"""
    MSplineBasis(x::AbstractVector{T};
                 boundary_knots::Union{Tuple{T,T},Nothing}=nothing,
                 interior_knots::Union{AbstractVector{T},Nothing}=nothing,
                 order::Int=4,
                 intercept::Bool=false,
                 df::Int=order - 1 + Int(intercept),
                 knots::Union{AbstractVector{T},Nothing}=nothing) where {T<:Real}

Calculate a basis for M-splines and return a callable type
for evaluating points in that basis.

# Keyword Arguments
- `boundary_knots`: boundary knots
- `interior_knots`: interior knots
- `order`: order of the spline
- `intercept`: bool for whether to include an intercept (column of ones). This behaviour is different to the Splines2 package from R, where intercept=FALSE will drop the first spline term.
- `df`: degrees of freedom
- `knots`: full set of knots (excluding repeats)

# Keyword Arguments of Returned Callable
- `center`: value to center the splines
- `derivs`: derivatives of the splines

# Examples

```jldoctest
julia> spline = MSplineBasis(collect(0.0:0.2:1.0); df=3);
julia> spline(collect(0.0:0.2:1.0))
6×3 Array{Float64,2}:
 0.0    0.0    0.0
 1.536  0.384  0.032
 1.728  1.152  0.256
 1.152  1.728  0.864
 0.384  1.536  2.048
 0.0    0.0    4.0
```
"""
function MSplineBasis(x::AbstractVector{T};
                      boundary_knots::Union{Tuple{T,T},Nothing}=nothing,
                      interior_knots::Union{AbstractVector{T},Nothing}=nothing,
                      order::Int=4,
                      intercept::Bool=false,
                      df::Int=order - 1 + Int(intercept),
                      knots::Union{AbstractVector{T},Nothing}=nothing) where {T<:Real}
    boundary_knots, interior_knots = spline_args(x,
                                                 boundary_knots,
                                                 interior_knots;
                                                 order,
                                                 intercept,
                                                 df,
                                                 knots)
    spline = MSplineBasis(boundary_knots, interior_knots, order, intercept)
    return spline
end

function (spline::MSplineBasis{T})(x::AbstractVector;
                                   derivs::Int=0,
                                   center::Union{Number,Nothing}=nothing) where {T}
    b = basis(spline, x, derivs)
    if !isnothing(center) && iszero(derivs)
        bc = basis(spline, T(center), derivs)
        for i in axes(b, 1)
            b[i, :] -= bc
        end
    end
    return b
end

"""
   ms(x::AbstractVector{T};
           center::Union{T,Nothing}=nothing,
           derivs::Int=0,
           kwargs...) where {T <: Real}

Calculate a basis for M-splines and return `x` expressed in that basis.

# Keyword arguments
- `center`: value to center the splines
- `derivs`: derivatives of the splines
- Further keyword arguments are passed to `MSplineBasis`.

# Arguments
- `boundary_knots :: Union{Tuple{T,T},Nothing} = nothing`: boundary knots
- `interior_knots :: Union{Array{T,1},Nothing} = nothing`: interior knots
- `order :: Int = 4`: order of the spline
- `intercept :: Bool = false`: bool for whether to include an intercept
- `df :: Int = order - 1 + Int(intercept)`: degrees of freedom
- `knots :: Union{Array{T,1}, Nothing} = nothing`: full set of knots (excluding repeats)

# Examples
```jldoctest
julia> ms(collect(0.0:0.2:1.0), df=3)
6×3 Array{Float64,2}:
 0.0    0.0    0.0
 1.536  0.384  0.032
 1.728  1.152  0.256
 1.152  1.728  0.864
 0.384  1.536  2.048
 0.0    0.0    4.0
```
"""
function ms(x::AbstractVector{T};
            center::Union{T,Nothing}=nothing,
            derivs::Int=0,
            kwargs...) where {T<:Real}
    spline = MSplineBasis(x; kwargs...)
    return spline(x; center, derivs)
end
