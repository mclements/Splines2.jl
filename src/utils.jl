zeroIndexedArray(a::Array{T}) where {T<:Real} = OffsetArray(a, Tuple(fill(-1, ndims(a))))

function _find_interval(bs::SplineBasis{T}, x::T) where {T<:Real}
    k = bs.order - 1
    n = bs.nknots - k - 1
    l = k
    while x < bs.knots[l] && l != k
        l = l - 1
    end
    l = l + 1
    while x >= bs.knots[l] && l != n
        l = l + 1
    end
    return l - 1
end

"""
    spline_args(x::AbstractVector{T};
                boundary_knots::Union{Tuple{T,T},Nothing}=nothing,
                interior_knots::Union{Vector{T},Nothing}=nothing,
                order::Int=4,
                intercept::Bool=false,
                df::Int=3 + Int(intercept),
                knots::Union{Vector{T},Nothing}=nothing,
                knots_offset::Int=0) where {T<:Real}

Utility function for processing the spline arguments.

Returns the computed boundary and interior knots. If these 
are already both computed (i.e. not nothing), then returns
the passed values.
"""
function spline_args(x::AbstractVector{T},
                     boundary_knots::Tuple{T,T},
                     interior_knots::Vector{T};
                     kwargs...) where {T<:Real}
    return (; boundary_knots, interior_knots)
end

function spline_args(x::AbstractVector{T},
                     boundary_knots::Union{Tuple{T,T},Nothing},
                     interior_knots::Union{Vector{T},Nothing};
                     order::Int=4,
                     intercept::Bool=false,
                     df::Int=3 + Int(intercept),
                     knots::Union{Vector{T},Nothing}=nothing,
                     knots_offset::Int=0) where {T<:Real}
    if !isnothing(knots)
        length(knots) == 1 &&
            error("At least two knots are required if knots are specified.")
        # if knots are passed, then use the first and last knots as
        # the boundary knots and any remaining knots as interior knots
        # NOTE: this assumes that the knots are already sorted so that the boundary knots
        # are first and last
        boundary_knots = (first(knots), last(knots))
        interior_knots = length(knots) == 2 ? nothing : knots[2:(length(knots) - 1)]
    else
        boundary_knots = @something(boundary_knots, extrema(x))
        iKnots = df - order + knots_offset + 1 - Int(intercept)
        if iKnots > 0
            # we exclude the 0th and 100th percentiles
            p = view(range(T(0); length=iKnots + 2, stop=T(1)), 2:(iKnots + 1))
            interior = boundary_knots[1] .<= x .<= boundary_knots[2]
            interior_knots = quantile(view(x, interior), p)
        end
    end
    return (; boundary_knots, interior_knots)
end
