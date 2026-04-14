struct ISplineBasis{T<:Real} <: AbstractSplineBasis{T}
    b_spline_basis::BSplineBasis{T}  # built with order+1, intercept=false
    order::Int                        # original order (b_spline_basis uses order+1)
    intercept::Bool
    interior_knots::Union{Vector{T},Nothing}
end

function ISplineBasis(boundary_knots::Tuple{T,T},
                      interior_knots::Union{Array{T,1},Nothing}=nothing,
                      order::Int=4,
                      intercept::Bool=false) where {T<:Real}
    spline = BSplineBasis(boundary_knots, interior_knots, order + 1, false)
    return ISplineBasis(spline, order, intercept, interior_knots)
end

function basis(is::ISplineBasis{T}, x::T, derivs::Int=0) where {T<:Real}
    vec = basis(is.b_spline_basis, x, derivs)
    order = is.order
    knots = parent(is.b_spline_basis.spline_basis.knots)
    js = isnothing(is.interior_knots) ? (order + 1) : searchsortedlast(knots, x)
    ncol = length(vec)
    for j in ncol:-1:1
        if j > js
            vec[j] = T(0)
        elseif j < ncol
            vec[j] = vec[j] + vec[j + 1]
        end
    end
    for j in ncol:-1:1
        if j < js - (order + 1) && iszero(derivs)
            vec[j] = T(1)
        end
    end
    return vec
end

"""
    ISplineBasis(x::AbstractVector{T};
                 boundary_knots::Union{Tuple{T,T},Nothing}=nothing,
                 interior_knots::Union{AbstractVector{T},Nothing}=nothing,
                 order::Int=4,
                 intercept::Bool=false,
                 df::Int=order + Int(intercept),
                 knots::Union{AbstractVector{T},Nothing}=nothing) where {T<:Real}

Calculate a basis for I-splines and return a callable type
for evaluating points in that basis.

# Keyword Arguments
- `boundary_knots`: boundary knots
- `interior_knots`: interior knots
- `order`: order of the spline
- `intercept`: bool for whether to include an intercept (column of ones). This behaviour is different to the splines2 package from R, where intercept=FALSE will drop the first spline term.
- `df`: degrees of freedom
- `knots`: full set of knots (excluding repeats)

# Keyword Arguments of Returned Callable
- `derivs`: derivatives of the splines

# Examples

```jldoctest
julia> spline = ISplineBasis(collect(0.0:0.2:1.0); df=3);
julia> spline(collect(0.0:0.2:1.0))
6×3 Array{Float64,2}:
 0.0     0.0     0.0
 0.1808  0.0272  0.0016
 0.5248  0.1792  0.0256
 0.8208  0.4752  0.1296
 0.9728  0.8192  0.4096
 1.0     1.0     1.0
```
"""
function ISplineBasis(x::AbstractVector{T};
                      boundary_knots::Union{Tuple{T,T},Nothing}=nothing,
                      interior_knots::Union{AbstractVector{T},Nothing}=nothing,
                      order::Int=4,
                      intercept::Bool=false,
                      df::Int=order + Int(intercept),
                      knots::Union{AbstractVector{T},Nothing}=nothing) where {T<:Real}
    boundary_knots, interior_knots = spline_args(x,
                                                 boundary_knots,
                                                 interior_knots;
                                                 order=order,
                                                 intercept,
                                                 df,
                                                 knots)
    spline = ISplineBasis(boundary_knots, interior_knots, order, intercept)
    return spline
end

function (spline::ISplineBasis{T})(x::AbstractVector; derivs::Int=0) where {T}
    b = basis(spline, x, derivs)
    if spline.intercept
        b = [ones(T, size(b, 1)) b]
    end
    return b
end

"""
   is(x::AbstractVector{T};
           derivs::Int=0,
           kwargs...) where {T <: Real}

Calculate a basis for I-splines and return `x` expressed in that basis.

# Keyword arguments
- `derivs`: derivatives of the splines
- Further keyword arguments are passed to `ISplineBasis`.

# Arguments
- `boundary_knots :: Union{Tuple{T,T},Nothing} = nothing`: boundary knots
- `interior_knots :: Union{Array{T,1},Nothing} = nothing`: interior knots
- `order :: Int = 4`: order of the spline
- `intercept :: Bool = false`: bool for whether to include an intercept
- `df :: Int = order + Int(intercept)`: degrees of freedom
- `knots :: Union{Array{T,1}, Nothing} = nothing`: full set of knots (excluding repeats)

# Examples
```jldoctest
julia> is(collect(0.0:0.2:1.0), df=3)
6×3 Array{Float64,2}:
 0.0     0.0     0.0
 0.1808  0.0272  0.0016
 0.5248  0.1792  0.0256
 0.8208  0.4752  0.1296
 0.9728  0.8192  0.4096
 1.0     1.0     1.0
```
"""
function is(x::AbstractVector{T}; derivs::Int=0, kwargs...) where {T<:Real}
    spline = ISplineBasis(x; kwargs...)
    return spline(x; derivs)
end
