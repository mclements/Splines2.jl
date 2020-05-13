module Splines2

export ns, ns_, bs, bs_, is, is_, ms, ms_, basis

using OffsetArrays
using LinearAlgebra
using Statistics

# utility functions
zeroArray(a :: Array{T}) where T<:Real = OffsetArray(a,(size(a).*0).-1)

# type definitions
abstract type AbstractSplineBasis{T<:Real} end
# Note: we have used OffsetArray for converting from C code
mutable struct SplineBasis{T<:Real} <: AbstractSplineBasis{T}
    order::Int # order of the spline
    nknots::Int # number of knots
    ncoef::Int # number of coefficients
    knots::OffsetArray{T,1} # knot vector
end
mutable struct BSplineBasis{T<:Real} <: AbstractSplineBasis{T}
    spline_basis::SplineBasis{T}
    boundary_knots::Tuple{T,T}
    interior_knots::Union{Array{T,1}, Nothing}
    intercept::Bool
    df::Int
end
mutable struct NSplineBasis{T<:Real} <: AbstractSplineBasis{T}
    b_spline_basis::BSplineBasis{T}
    qmat::Array{T,2}
    tl0::Array{T,1}
    tl1::Array{T,1}
    tr0::Array{T,1}
    tr1::Array{T,1}
end

# constructors
function SplineBasis(knots :: Array{T,1}, order :: Int = 4) where T<:Real
    SplineBasis(order,
                length(knots),
                length(knots)-order, 
                zeroArray(knots))
end
function BSplineBasis(boundary_knots :: Tuple{T,T},
                      interior_knots :: Union{Array{T,1}, Nothing} = nothing,
                      order :: Int = 4,
                      intercept :: Bool = false) where T<:Real
    l_interior_knots = interior_knots == nothing ? 0 : length(interior_knots)
    df = Int(intercept) + order - 1 + l_interior_knots
    nknots = l_interior_knots + 2*order
    ncoef = nknots - order
    knots = zeros(T,nknots)
    for i=1:order
        knots[i] = boundary_knots[1]
        knots[nknots-i+1] = boundary_knots[2]
    end
    if (l_interior_knots>0)
        for i=1:l_interior_knots
            knots[i+order]=interior_knots[i]
        end
    end
    BSplineBasis(SplineBasis(knots,order), boundary_knots, interior_knots, intercept, df)
end
function NSplineBasis(boundary_knots :: Tuple{T,T},
                      interior_knots :: Union{Array{T,1}, Nothing} = nothing,
                      order :: Int = 4,
                      intercept :: Bool = false) :: NSplineBasis{T} where T<:Real
    bs = BSplineBasis(boundary_knots,interior_knots,order,intercept)
    co = basis(bs,[bs.boundary_knots[1],bs.boundary_knots[2]],2)
    qmat_ = LinearAlgebra.qr(transpose(co)).Q
    qmat = collect(transpose(qmat_[:, 3:size(qmat_,2)]))
    tl0 = qmat * basis(bs,boundary_knots[1])
    tl1 = qmat * basis(bs,boundary_knots[1],1)
    tr0 = qmat * basis(bs,boundary_knots[2])
    tr1 = qmat * basis(bs,boundary_knots[2],1)
    NSplineBasis(bs,
                 qmat,
                 tl0,tl1,tr0,tr1)
end

function basis(bs :: SplineBasis{T}, x :: T, ders :: Int = 0) where T<:Real
    function find_interval(x :: T)
        k = bs.order - 1
        n = bs.nknots - k - 1
        l = k
        while (x < bs.knots[l] && l != k)
            l = l - 1
        end
        l = l + 1
        while (x >= bs.knots[l] && l != n)
            l = l + 1
        end
        l - 1
    end
    t = bs.knots
    k = bs.order - 1
    m = ders
    hh = bs.order
    ell = find_interval(x)
    result = zeroArray(zeros(T,2*k+2))
    one = T(1)
    zero = T(0)
    result[0] = one
    for j = 1:(k-m)
        for n = 0:(j-1)
            result[hh+n] = result[n]
        end
        result[0] = zero
        for n = 1:j
            ind = ell + n
            xb = t[ind]
            xa = t[ind - j]
            if (xb == xa) 
                result[n] = zero
                continue
            end
            w = result[hh + n - 1]/(xb - xa)
            result[n - 1] = result[n-1] + w*(xb - x)
            result[n] = w*(x - xa)
        end
    end
    for j = (k - m + 1):k
        for n = 0:(j-1)
            result[hh+n] = result[n]
        end
        result[0] = zero
        for n = 1:j
            ind = ell + n
            xb = t[ind]
            xa = t[ind - j]
            if (xb == xa)
                result[m] = zero
                continue
            end
            w = j*result[hh + n - 1]/(xb - xa)
            result[n - 1] = result[n - 1] - w
            result[n] = w
        end
    end
    offset = ell - k
    v = zeros(T,bs.ncoef)
    v[(1+offset):(k+1+offset)] = result[0:k]
    v
end

function basis(bs :: BSplineBasis{T}, x :: T, ders :: Int = 0) where T<:Real
    if (x < bs.boundary_knots[1] || x > bs.boundary_knots[2])
        if (x < bs.boundary_knots[1])
            k_pivot = T(0.75)*bs.boundary_knots[1]+
            T(0.25)*bs.spline_basis.knots[bs.spline_basis.order+1-1] # 0-based
        else
            k_pivot = T(0.75)*bs.boundary_knots[2]+
            T(0.25)*bs.spline_basis.knots[length(bs.spline_basis.knots)-bs.spline_basis.order-1-1] # 0-based
        end
        delta = x - k_pivot
        if (ders==0)
            vec = basis(bs.spline_basis,k_pivot,0) +
                basis(bs.spline_basis,k_pivot,1)*delta +
                basis(bs.spline_basis,k_pivot,2)*delta*delta/T(2.0) + 
                basis(bs.spline_basis,k_pivot,3)*delta*delta*delta/T(6.0)
        elseif (ders==1)
            vec = splines.basis(bs.spline_basis,k_pivot,1) +
                splines.basis(bs.spline_basis,k_pivot,2)*delta + 
                splines.basis(bs.spline_basis,k_pivot,3)*delta*delta/T(2.0)
        elseif (ders==2)
            vec = splines.basis(bs.spline_basis,k_pivot,2) + 
                splines.basis(bs.spline_basis,k_pivot,3)*delta
        elseif (ders==3)
            vec = splines.basis(bs.spline_basis,k_pivot,3)
        else
            vec = k_pivot .* T(0)
        end
    else 
        vec = basis(bs.spline_basis, x, ders)
    end
    if (!bs.intercept)
        vec = vec[2:length(vec)]
    end
    vec
end

function basis(ns :: NSplineBasis{T}, x::T, ders :: Int = 0) where T<:Real
    if (x < ns.b_spline_basis.boundary_knots[1])
        if (ders==0)
            vec = ns.tl0 + (x - ns.b_spline_basis.boundary_knots[1])*ns.tl1
        elseif (ders==1)
            vec = ns.tl1
        else
            vec = ns.tl1.* T(0)
        end
    elseif (x > ns.b_spline_basis.boundary_knots[2])
        if (ders==0)
            vec = ns.tr0 + (x - ns.b_spline_basis.boundary_knots[2])*ns.tr1
        elseif (ders==1)
            vsc = ns.tr1
        else
            vec = ns.tr1.* T(0)
        end
    else 
        vec = ns.qmat*basis(ns.b_spline_basis, x, ders)
    end
    vec
end

function basis(bs :: AbstractSplineBasis{T}, x :: Array{T,1}, ders :: Int = 0) where T<:Real
    f(xi) = basis(bs, xi, ders)
    copy(transpose(reduce(hcat, f.(x))))
end

# utility function for processing the spline arguments
function spline_args(x :: Array{T,1};
                     boundary_knots :: Union{Tuple{T,T},Nothing} = nothing,
                     interior_knots :: Union{Array{T,1},Nothing} = nothing,
                     order :: Int = 4,
                     intercept :: Bool = false,
                     df :: Int = 3 + Int(intercept),
                     knots :: Union{Array{T,1}, Nothing} = nothing,
                     knots_offset :: Int = 0) where T<:Real
    if (interior_knots != nothing && boundary_knots != nothing)
        # pass
    elseif (knots != nothing)
        boundary_knots = extrema(knots)
        interior_knots = length(knots)==2 ? nothing : knots[2:(length(knots)-1)]
    else
        if (boundary_knots == nothing)
            boundary_knots = extrema(x)
        end
        iKnots = df - order + knots_offset + 1 - Int(intercept)
        if (iKnots>0)
            p = range(T(0), length=iKnots+2, stop=T(1))[2:(iKnots+1)]
            index = (x .>= boundary_knots[1]) .* (x .<= boundary_knots[2])
            interior_knots = Statistics.quantile(x[index], p)
        end
    end
    (boundary_knots, interior_knots)
end

"""
    bs_(x :: Array{T,1}; <keyword arguments>) where T<:Real

Calculate a basis for B-splines and return a function with signature
`(x:: Array{T,1}; ders :: Int = 0)` for evaluation of `ders`
derivative for the splines at `x`.

The keyword arguments include one of:
1. `df`, possibly in combination with `intercept`
2. `boundary_knots` and `interior_knots`
3. `knots`

# Arguments
- `boundary_knots :: Union{Tuple{T,T},Nothing} = nothing`: boundary knots
- `interior_knots :: Union{Array{T,1},Nothing} = nothing`: interior knots
- `order :: Int = 4`: order of the spline
- `intercept :: Bool = false`: bool for whether to include an intercept
- `df :: Int = order - 1 + Int(intercept)`: degrees of freedom
- `knots :: Union{Array{T,1}, Nothing} = nothing`: full set of knots (excluding repeats)
- `centre :: Union{T,Nothing} = nothing)`: value to centre the splines

# Examples
```jldoctest
julia> Splines2.bs_(collect(0.0:0.2:1.0), df=3)(collect(0.0:0.2:1.0))
6×3 Array{Float64,2}:
 0.0    0.0    0.0  
 0.384  0.096  0.008
 0.432  0.288  0.064
 0.288  0.432  0.216
 0.096  0.384  0.512
 0.0    0.0    1.0  
```
"""
function bs_(x :: Array{T,1};
            boundary_knots :: Union{Tuple{T,T},Nothing} = nothing,
            interior_knots :: Union{Array{T,1},Nothing} = nothing,
            order :: Int = 4,
            intercept :: Bool = false,
            df :: Int = order - 1 + Int(intercept),
            knots :: Union{Array{T,1}, Nothing} = nothing,
            centre :: Union{T,Nothing} = nothing) where T<:Real
    (boundary_knots, interior_knots) =
        spline_args(x, boundary_knots=boundary_knots, interior_knots=interior_knots,
                    order=order, intercept=intercept,
                    df=df, knots=knots)
    spline = BSplineBasis(boundary_knots, interior_knots, order, intercept)
    function eval(x :: Array{T,1}; ders :: Int = 0)
        b = basis(spline, x, ders)
        if (centre != nothing && ders==0)
            bc = basis(spline, centre, ders)
            for i=1:size(b,1)
                b[i,:] -= bc
            end
        end
        b
    end
    eval
end

"""
    bs(x :: Array{T,1}; <keyword arguments>) where T<:Real

Calculate a basis for B-splines. 

The keyword arguments include one of:
1. `df`, possibly in combination with `intercept`
2. `boundary_knots` and `interior_knots`
3. `knots`

# Arguments
- `boundary_knots :: Union{Tuple{T,T},Nothing} = nothing`: boundary knots
- `interior_knots :: Union{Array{T,1},Nothing} = nothing`: interior knots
- `order :: Int = 4`: order of the spline
- `intercept :: Bool = false`: bool for whether to include an intercept
- `df :: Int = order - 1 + Int(intercept)`: degrees of freedom
- `knots :: Union{Array{T,1}, Nothing} = nothing`: full set of knots (excluding repeats)
- `centre :: Union{T,Nothing} = nothing)`: value to centre the splines
- `ders :: Int = 0`: derivatives of the splines

# Examples
```jldoctest
julia> Splines2.bs(collect(0.0:0.2:1.0), df=3)
6×3 Array{Float64,2}:
 0.0    0.0    0.0  
 0.384  0.096  0.008
 0.432  0.288  0.064
 0.288  0.432  0.216
 0.096  0.384  0.512
 0.0    0.0    1.0  
```
"""
function bs(x :: Array{T,1}; ders :: Int = 0, kwargs...) where T<:Real
    bs_(x; kwargs...)(x, ders=ders)
end

"""
    ns_(x :: Array{T,1}; <keyword arguments>) where T<:Real

Calculate a basis for natural B-splines and return a function with signature
`(x:: Array{T,1}; ders :: Int = 0)` for evaluation of `ders`
derivative for the splines at `x`.

The keyword arguments include one of:
1. `df`, possibly in combination with `intercept`
2. `boundary_knots` and `interior_knots`
3. `knots`

# Arguments
- `boundary_knots :: Union{Tuple{T,T},Nothing} = nothing`: boundary knots
- `interior_knots :: Union{Array{T,1},Nothing} = nothing`: interior knots
- `order :: Int = 4`: order of the spline
- `intercept :: Bool = false`: bool for whether to include an intercept
- `df :: Int = order - 1 + Int(intercept)`: degrees of freedom
- `knots :: Union{Array{T,1}, Nothing} = nothing`: full set of knots (excluding repeats)
- `centre :: Union{T,Nothing} = nothing)`: value to centre the splines

# Examples
```jldoctest
julia> Splines2.ns_(collect(0.0:0.2:1.0), df=3)(collect(0.0:0.2:1.0))
6×3 Array{Float64,2}:
  0.0       0.0        0.0     
 -0.100444  0.409332  -0.272888
  0.102383  0.540852  -0.359235
  0.501759  0.386722  -0.172481
  0.418872  0.327383   0.217745
 -0.142857  0.428571   0.714286
```
"""
function ns_(x :: Array{T,1};
            boundary_knots :: Union{Tuple{T,T},Nothing} = nothing,
            interior_knots :: Union{Array{T,1},Nothing} = nothing,
            order :: Int = 4,
            intercept :: Bool = false,
            df :: Int = order - 3 + Int(intercept),
            knots :: Union{Array{T,1}, Nothing} = nothing,
            centre :: Union{T,Nothing} = nothing) where T<:Real
    (boundary_knots, interior_knots) =
        spline_args(x, boundary_knots=boundary_knots, interior_knots=interior_knots,
                    order=order, intercept=intercept, df=df, knots=knots, knots_offset=2)
    spline = NSplineBasis(boundary_knots, interior_knots, order, intercept)
    function eval(x :: Array{T,1}; ders :: Int = 0)
        b = basis(spline, x, ders)
        if (centre != nothing && ders==0)
            bc = basis(spline, centre, ders)
            for i=1:size(b,1)
                b[i,:] -= bc
            end
        end
        b
    end
    eval
end

"""
    ns(x :: Array{T,1}; <keyword arguments>) where T<:Real

Calculate a basis for natural B-splines. 

The keyword arguments include one of:
1. `df`, possibly in combination with `intercept`
2. `boundary_knots` and `interior_knots`
3. `knots`

# Arguments
- `boundary_knots :: Union{Tuple{T,T},Nothing} = nothing`: boundary knots
- `interior_knots :: Union{Array{T,1},Nothing} = nothing`: interior knots
- `order :: Int = 4`: order of the spline
- `intercept :: Bool = false`: bool for whether to include an intercept
- `df :: Int = order - 1 + Int(intercept)`: degrees of freedom
- `knots :: Union{Array{T,1}, Nothing} = nothing`: full set of knots (excluding repeats)
- `centre :: Union{T,Nothing} = nothing)`: value to centre the splines
- `ders :: Int = 0`: derivatives of the splines

# Examples
```jldoctest
julia> Splines2.ns(collect(0.0:0.2:1.0), df=3)
6×3 Array{Float64,2}:
  0.0       0.0        0.0     
 -0.100444  0.409332  -0.272888
  0.102383  0.540852  -0.359235
  0.501759  0.386722  -0.172481
  0.418872  0.327383   0.217745
 -0.142857  0.428571   0.714286
```
"""
function ns(x :: Array{T,1}; ders :: Int = 0, kwargs...) where T<:Real
    ns_(x; kwargs...)(x, ders=ders)
end

"""
    is_(x :: Array{T,1}; <keyword arguments>) where T<:Real

Calculate a basis for I-splines and return a function with signature
`(x:: Array{T,1}; ders :: Int = 0)` for evaluation of `ders`
derivative for the splines at `x`.

The keyword arguments include one of:
1. `df`, possibly in combination with `intercept`
2. `boundary_knots` and `interior_knots`
3. `knots`

# Arguments
- `boundary_knots :: Union{Tuple{T,T},Nothing} = nothing`: boundary knots
- `interior_knots :: Union{Array{T,1},Nothing} = nothing`: interior knots
- `order :: Int = 4`: order of the spline
- `intercept :: Bool = false`: bool for whether to include an intercept
- `df :: Int = order - 1 + Int(intercept)`: degrees of freedom
- `knots :: Union{Array{T,1}, Nothing} = nothing`: full set of knots (excluding repeats)

# Examples
```jldoctest
julia> Splines2.is_(collect(0.0:0.2:1.0), df=3)(collect(0.0:0.2:1.0))
6×3 Array{Float64,2}:
 0.0     0.0     0.0   
 0.1808  0.0272  0.0016
 0.5248  0.1792  0.0256
 0.8208  0.4752  0.1296
 0.9728  0.8192  0.4096
 1.0     1.0     1.0   
```
"""
function is_(x :: Array{T,1};
            boundary_knots :: Union{Tuple{T,T},Nothing} = nothing,
            interior_knots :: Union{Array{T,1},Nothing} = nothing,
            order :: Int = 4,
            intercept :: Bool = false,
            df :: Int = order - 1 + Int(intercept),
            knots :: Union{Array{T,1}, Nothing} = nothing) where T<:Real
    (boundary_knots, interior_knots) =
        spline_args(x, boundary_knots=boundary_knots, order=order+1,
                    intercept=intercept, df=df+1, knots=knots)
    spline = BSplineBasis(boundary_knots, interior_knots, order+1, false)
    knots = parent(spline.spline_basis.knots)
    findj(x) = interior_knots==nothing ? (order+1) : searchsortedlast(knots,x)
    function eval(x :: Array{T,1}; ders :: Int = 0)
        b = basis(spline,x,ders)
        (nrow,ncol) = size(b)
        for i=1:nrow
            js = findj(x[i])
            for j in ncol:-1:1
                if j>js
                    b[i,j]=T(0)
                elseif j==ncol
                    b[i,j]=b[i,j]
                else
                    b[i,j]=b[i,j]+b[i,j+1]
                end
            end
            for j in ncol:-1:1
                if j < js-(order+1) && der==0
                    b[i,j]=T(1)
                end
            end
        end
        intercept ? b : b[:, 2:size(b,2)]
    end
    eval
end

"""
    is(x :: Array{T,1}; <keyword arguments>) where T<:Real

Calculate a basis for I-splines. 

The keyword arguments include one of:
1. `df`, possibly in combination with `intercept`
2. `boundary_knots` and `interior_knots`
3. `knots`

# Arguments
- `boundary_knots :: Union{Tuple{T,T},Nothing} = nothing`: boundary knots
- `interior_knots :: Union{Array{T,1},Nothing} = nothing`: interior knots
- `order :: Int = 4`: order of the spline
- `intercept :: Bool = false`: bool for whether to include an intercept
- `df :: Int = order - 1 + Int(intercept)`: degrees of freedom
- `knots :: Union{Array{T,1}, Nothing} = nothing`: full set of knots (excluding repeats)
- `ders :: Int = 0`: derivatives of the splines

# Examples
```jldoctest
julia> Splines2.is(collect(0.0:0.2:1.0), df=3)
6×3 Array{Float64,2}:
 0.0     0.0     0.0   
 0.1808  0.0272  0.0016
 0.5248  0.1792  0.0256
 0.8208  0.4752  0.1296
 0.9728  0.8192  0.4096
 1.0     1.0     1.0   
```
"""
function is(x :: Array{T,1}; ders :: Int = 0, kwargs...) where T<:Real
    is_(x; kwargs...)(x, ders=ders)
end

"""
    ms_(x :: Array{T,1}; <keyword arguments>) where T<:Real

Calculate a basis for M-splines and return a function with signature
`(x:: Array{T,1}; ders :: Int = 0)` for evaluation of `ders`
derivative for the splines at `x`.

The keyword arguments include one of:
1. `df`, possibly in combination with `intercept`
2. `boundary_knots` and `interior_knots`
3. `knots`

# Arguments
- `boundary_knots :: Union{Tuple{T,T},Nothing} = nothing`: boundary knots
- `interior_knots :: Union{Array{T,1},Nothing} = nothing`: interior knots
- `order :: Int = 4`: order of the spline
- `intercept :: Bool = false`: bool for whether to include an intercept
- `df :: Int = order - 1 + Int(intercept)`: degrees of freedom
- `knots :: Union{Array{T,1}, Nothing} = nothing`: full set of knots (excluding repeats)
- `centre :: Union{T,Nothing} = nothing)`: value to centre the splines

# Examples
```jldoctest
julia> Splines2.ms_(collect(0.0:0.2:1.0), df=3)(collect(0.0:0.2:1.0))
6×3 Array{Float64,2}:
 0.0    0.0    0.0  
 1.536  0.384  0.032
 1.728  1.152  0.256
 1.152  1.728  0.864
 0.384  1.536  2.048
 0.0    0.0    4.0  
```
"""
function ms_(x :: Array{T,1};
            boundary_knots :: Union{Tuple{T,T},Nothing} = nothing,
            interior_knots :: Union{Array{T,1},Nothing} = nothing,
            order :: Int = 4,
            intercept :: Bool = false,
            df :: Int = order - 1 + Int(intercept),
            knots :: Union{Array{T,1}, Nothing} = nothing,
            centre :: Union{T,Nothing} = nothing) where T<:Real
    (boundary_knots, interior_knots) =
        spline_args(x, boundary_knots=boundary_knots, interior_knots=interior_knots,
                    order=order, intercept=intercept, df=df, knots=knots)
    spline = BSplineBasis(boundary_knots, interior_knots, order, intercept)
    function eval(x :: Array{T,1}; ders :: Int = 0)
        b = basis(spline, x, ders)
        knots = parent(spline.spline_basis.knots)
        function loop(j)
            denom = knots[j+order]-knots[j]
            denom > T(0) ? order/denom : T(0)
        end
        transCoef = map(loop,1:(length(knots)-order))
        if (!intercept)
            transCoef = transCoef[2:length(transCoef)]
        end
        if (centre != nothing && ders==0)
            bc = basis(spline, centre, ders)
            for i=1:size(b,1)
                b[i,:] -= bc
            end
        end
        for i=1:size(b,1)
            b[i,:] .*= transCoef 
        end
        b
    end
    eval
end

"""
    ms(x :: Array{T,1}; <keyword arguments>) where T<:Real

Calculate a basis for M-splines. 

The keyword arguments include one of:
1. `df`, possibly in combination with `intercept`
2. `boundary_knots` and `interior_knots`
3. `knots`

# Arguments
- `boundary_knots :: Union{Tuple{T,T},Nothing} = nothing`: boundary knots
- `interior_knots :: Union{Array{T,1},Nothing} = nothing`: interior knots
- `order :: Int = 4`: order of the spline
- `intercept :: Bool = false`: bool for whether to include an intercept
- `df :: Int = order - 1 + Int(intercept)`: degrees of freedom
- `knots :: Union{Array{T,1}, Nothing} = nothing`: full set of knots (excluding repeats)
- `centre :: Union{T,Nothing} = nothing)`: value to centre the splines
- `ders :: Int = 0`: derivatives of the splines

# Examples
```jldoctest
julia> Splines2.ms(collect(0.0:0.2:1.0), df=3)
6×3 Array{Float64,2}:
 0.0    0.0    0.0  
 1.536  0.384  0.032
 1.728  1.152  0.256
 1.152  1.728  0.864
 0.384  1.536  2.048
 0.0    0.0    4.0  
```
"""
function ms(x :: Array{T,1}; ders :: Int = 0, kwargs...) where T<:Real
    ms_(x; kwargs...)(x, ders=ders)
end

end # module
