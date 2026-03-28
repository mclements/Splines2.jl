struct NSplineBasis{T<:Real} <: AbstractSplineBasis{T}
    b_spline_basis::BSplineBasis{T}
    qmat::Matrix{T}
    tl0::Vector{T}
    tl1::Vector{T}
    tr0::Vector{T}
    tr1::Vector{T}
end

function NSplineBasis(boundary_knots::Tuple{T,T},
                      interior_knots::Union{Array{T,1},Nothing}=nothing,
                      order::Int=4,
                      intercept::Bool=false)::NSplineBasis{T} where {T<:Real}
    bs = BSplineBasis(boundary_knots, interior_knots, order, intercept)
    co = basis(bs, [bs.boundary_knots[1], bs.boundary_knots[2]], 2)
    qmat_ = LinearAlgebra.qr(transpose(co)).Q
    qmat = Matrix(transpose(@view qmat_[:, 3:size(qmat_, 2)]))
    tl0 = qmat * basis(bs, boundary_knots[1])
    tl1 = qmat * basis(bs, boundary_knots[1], 1)
    tr0 = qmat * basis(bs, boundary_knots[2])
    tr1 = qmat * basis(bs, boundary_knots[2], 1)
    return NSplineBasis(bs, qmat, tl0, tl1, tr0, tr1)
end

function basis(ns::NSplineBasis{T}, x::T, ders::Int=0) where {T<:Real}
    if (x < ns.b_spline_basis.boundary_knots[1])
        if (ders == 0)
            vec = ns.tl0 + (x - ns.b_spline_basis.boundary_knots[1]) * ns.tl1
        elseif (ders == 1)
            vec = ns.tl1
        else
            vec = ns.tl1 .* T(0)
        end
    elseif (x > ns.b_spline_basis.boundary_knots[2])
        if (ders == 0)
            vec = ns.tr0 + (x - ns.b_spline_basis.boundary_knots[2]) * ns.tr1
        elseif (ders == 1)
            vec = ns.tr1
        else
            vec = ns.tr1 .* T(0)
        end
    else
        vec = ns.qmat * basis(ns.b_spline_basis, x, ders)
    end
    return vec
end

"""
    NSplineBasis(x::AbstractVector{T};              
                 boundary_knots::Union{Tuple{T,T},Nothing}=nothing,
                 interior_knots::Union{AbstractVector{T},Nothing}=nothing,
                 order::Int=4,
                 intercept::Bool=false,
                 df::Int=order - 3 + Int(intercept),
                 knots::Union{AbstractVector{T},Nothing}=nothing) where {T<:Real}

Calculate a basis for natural B-splines and return a callable type
for evaluating points in that basis.

# Keyword Arguments
- `boundary_knots`: boundary knots
- `interior_knots`: interior knots
- `order`: order of the spline
- `intercept`: bool for whether to include an intercept
- `df`: degrees of freedom
- `knots`: full set of knots (excluding repeats)

# Keyword Arguments of Returned Callable
- `center`: value to center the splines
- `derivs`: derivatives of the splines

# Examples

```jldoctest
julia> spline = NSplineBasis(x; df=3);
julia> spline(0.0:0.2:1.0)
6×3 Matrix{Float64}:
  0.0       0.0        0.0     
 -0.100444  0.409332  -0.272888
  0.102383  0.540852  -0.359235
  0.501759  0.386722  -0.172481
  0.418872  0.327383   0.217745
 -0.142857  0.428571   0.714286
```
"""
function NSplineBasis(x::AbstractVector{T};
                      boundary_knots::Union{Tuple{T,T},Nothing}=nothing,
                      interior_knots::Union{AbstractVector{T},Nothing}=nothing,
                      order::Int=4,
                      intercept::Bool=false,
                      df::Int=order - 3 + Int(intercept),
                      knots::Union{AbstractVector{T},Nothing}=nothing) where {T<:Real}
    boundary_knots, interior_knots = spline_args(x,
                                                 boundary_knots,
                                                 interior_knots;
                                                 order,
                                                 intercept,
                                                 df,
                                                 knots,
                                                 knots_offset=2)
    spline = NSplineBasis(boundary_knots, interior_knots, order, intercept)
    return spline
end

function (spline::NSplineBasis{T})(x::AbstractVector;
                                   derivs::Int=0,
                                   center::Union{Number,Nothing}=nothing) where {T}
    b = basis(spline, x, derivs)
    if !isnothing(center) && iszero(derivs)
        bc = basis(spline, center, derivs)
        for i in axes(b, 1)
            b[i, :] -= bc
        end
    end
    return b
end

"""
   ns(x::AbstractVector{T};            
            center::Union{T,Nothing}=nothing,
            derivs::Int=0,
            kwargs...) where {T <: Real}

Calculate a basis for natural B-splines and return `x` expressed in that basis. 

# Keyword arguments
- `center`: value to center the splines
- `derivs`: derivatives of the splines
- Further keyword arguments are passed to `NSplineBasis`.

# Arguments
- `boundary_knots :: Union{Tuple{T,T},Nothing} = nothing`: boundary knots
- `interior_knots :: Union{Array{T,1},Nothing} = nothing`: interior knots
- `order :: Int = 4`: order of the spline
- `intercept :: Bool = false`: bool for whether to include an intercept
- `df :: Int = order - 1 + Int(intercept)`: degrees of freedom
- `knots :: Union{Array{T,1}, Nothing} = nothing`: full set of knots (excluding repeats)

# Examples
```jldoctest
julia> ns(0.0:0.2:1.0; df=3)

6×3 Matrix{Float64}:
  0.0       0.0        0.0     
 -0.100444  0.409332  -0.272888
  0.102383  0.540852  -0.359235
  0.501759  0.386722  -0.172481
  0.418872  0.327383   0.217745
 -0.142857  0.428571   0.714286
```
"""
function ns(x::AbstractVector{T};
            center::Union{T,Nothing}=nothing,
            derivs::Int=0,
            kwargs...) where {T<:Real}
    spline = NSplineBasis(x; kwargs...)

    return spline(x; center, derivs)
end

# TODO add show methods
