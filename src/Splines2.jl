module Splines2

export ns, ns_, bs, bs_, is, is_, ms, ms_

using OffsetArrays
using LinearAlgebra
using Statistics

# utility functions
zeroArray(a :: Array{T}) where T<:Real = OffsetArray(a,(size(a).*0).-1)

# type definitions
abstract type AbstractSplineBasis{T<:Real} end
# Note: we have used OffsetArray for converting from C code
mutable struct SplineBasis{T<:Real} <: AbstractSplineBasis{T}
    order::Int32 # order of the spline
    ordm1::Int32 # order -1 (3 for cubic splines
    nknots::Int32 # number of knots
    curs::Int32 # current position in knots vector
    boundary::Int32 
    ncoef::Int32 # number of coefficients
    ldel::OffsetArray{T,1} # differences from knots on the left
    rdel::OffsetArray{T,1} # differences from knots on the right
    knots::OffsetArray{T,1} # knot vector
    a::OffsetArray{T,1} # scratch vector
end
mutable struct BSplineBasis{T<:Real} <: AbstractSplineBasis{T}
    spline_basis::SplineBasis{T}
    boundary_knots::Tuple{T,T}
    interior_knots::Union{Array{T,1}, Nothing}
    intercept::Bool
    df::Int32
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
function SplineBasis(knots :: Array{T,1}, order :: Int32 = 4) where T<:Real
    ordm1 = order -1 
    SplineBasis(order,
                ordm1,
                length(knots),
                0,
                0,
                length(knots)-order, 
                zeroArray(zeros(T,ordm1)),
                zeroArray(zeros(T,ordm1)),
                zeroArray(knots),
                zeroArray(zeros(T,order)))
end
function BSplineBasis(boundary_knots :: Tuple{T,T},
                      interior_knots :: Union{Array{T,1}, Nothing} = nothing,
                      order :: Int32 = 4,
                      intercept :: Bool = false) where T<:Real
    l_interior_knots = interior_knots == nothing ? 0 : length(interior_knots)
    df = Int32(intercept) + order - 1 + l_interior_knots
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
                      order :: Int32 = 4,
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

function basis(bs :: SplineBasis{T}, x :: T, ders :: Int32 = 0) where T<:Real
    function set_cursor!(x :: T)
        bs.curs = -1
        bs.boundary = T(0)
        for i = 0:(bs.nknots - 1)
	    if (bs.knots[i] >= x)
                bs.curs = i
            end
	    if (bs.knots[i] > x)
                break
            end
        end
        if (bs.curs > bs.nknots - bs.order) 
	    lastLegit = bs.nknots - bs.order
	    if (x == bs.knots[lastLegit]) 
	        bs.boundary = 1
                bs.curs = lastLegit
            end
        end
        return bs.curs
    end
    function diff_table!(x :: T, ndiff :: Int32)
        for i = 0:(ndiff-1)
	    bs.rdel[i] = bs.knots[bs.curs + i] - x
	    bs.ldel[i] = x - bs.knots[bs.curs - (i + 1)]
        end
    end
    function slow_evaluate(x :: T)
        ti    = bs.curs
        outer = bs.ordm1
        nder = ders # copy
        if (bs.boundary>0 && nder == bs.ordm1) # value is arbitrary
	    return T(0)
        end
        while(nder>0)   # FIXME: divides by zero
            nder -= 1
            apt = 0
            lpt = ti - outer
            for inner = (outer-1):-1:0 
	        bs.a[apt] = T(outer) * (bs.a[apt + 1] - bs.a[apt])/(bs.knots[lpt + outer] - bs.knots[lpt])
                apt += 1
                lpt += 1
            end
	    outer -= 1
        end
        diff_table!(x, outer)
        while(outer>0)
            outer -= 1
            apt = 0
            lpt = outer
            rpt = 0
            for inner = outer:-1:0
	        # FIXME: divides by zero
	        bs.a[apt] = (bs.a[apt + 1] * bs.ldel[lpt] + bs.a[apt] * bs.rdel[rpt])/(bs.rdel[rpt] + bs.ldel[lpt])
                lpt -= 1
                rpt += 1
                apt += 1
            end
        end
        return(bs.a[0])
    end
    # fast evaluation of basis functions
    function basis_funcs(x :: T)
        b = zeroArray(zeros(T,bs.order))
        diff_table!(x, bs.ordm1)
        b[0] = T(1)
        for j = 1:bs.ordm1
	    saved = T(0)
	    for r = 0:(j-1) # do not divide by zero
	        den = bs.rdel[r] + bs.ldel[j - 1 - r]
	        if(den != T(0)) 
	            term = b[r]/den
	            b[r] = saved + bs.rdel[r] * term
	            saved = bs.ldel[j - 1 - r] * term
	        else
	            if(r != T(0) || bs.rdel[r] != T(0))
	                b[r] = saved
                    end
	            saved = T(0)
	        end
	    end
	    b[j] = saved
        end
        return(b)
    end
    val = zeroArray(zeros(T,bs.ncoef))
    set_cursor!(x)
    io = bs.curs - bs.order
    if (io < 0 || io > bs.nknots) 
	for j=0:(bs.order-1)
	    val[j+io] = T(0) # R_NaN
	end
    elseif (ders > 0) # slow method for derivatives
	for i = 0:(bs.order-1)
	    for j=0:(bs.order-1)
                bs.a[j] = T(0)
            end
	    bs.a[i] = T(1)
	    val[i+io] = slow_evaluate(x)
	end
    else	# fast method for value
	valtmp = basis_funcs(x)
	for i=0:(length(valtmp)-1)
	    val[i+io]=valtmp[i]
        end
    end
    return parent(val)
end

function basis(bs :: BSplineBasis{T}, x :: T, ders :: Int32 = 0) where T<:Real
    if (x < bs.boundary_knots[1] || x > bs.boundary_knots[2])
        if (x < bs.boundary_knots[1])
            k_pivot = T(0.75)*bs.boundary_knots[1]+
            T(0.25)*bs.spline_basis.knots[5-1] # 0-based
        else
            k_pivot = T(0.75)*bs.boundary_knots[2]+
            T(0.25)*bs.spline_basis.knots[length(bs.spline_basis.knots)-4-1] # 0-based
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

function basis(ns :: NSplineBasis{T}, x::T, ders :: Int32 = 0) where T<:Real
    if (x < ns.b_spline_basis.boundary_knots[1])
        if (ders==0)
            return(ns.tl0 + (x - ns.b_spline_basis.boundary_knots[1])*ns.tl1)
        elseif (ders==1)
            return(ns.tl1)
        else
            return(ns.tl1.* T(0))
        end
    elseif (x > ns.b_spline_basis.boundary_knots[2])
        if (ders==0)
            return(ns.tr0 + (x - ns.b_spline_basis.boundary_knots[2])*ns.tr1)
        elseif (ders==1)
            return(ns.tr1)
        else
            return(ns.tr1.* T(0))
        end
    else 
        return(ns.qmat*basis(ns.b_spline_basis, x, ders))
    end
end

function basis(bs :: AbstractSplineBasis{T}, x :: Array{T,1}, ders :: Int32 = 0) where T<:Real
    f(xi) = basis(bs, xi, ders)
    copy(transpose(reduce(hcat, f.(x))))
end

# utility function for processing the spline arguments
function spline_args(x :: Array{T,1};
                     boundary_knots :: Union{Tuple{T,T},Nothing} = nothing,
                     interior_knots :: Union{Array{T,1},Nothing} = nothing,
                     order :: Int32 = 4,
                     intercept :: Bool = false,
                     df :: Int32 = 3 + Int32(intercept),
                     knots :: Union{Array{T,1}, Nothing} = nothing,
                     knots_offset :: Int32 = 0) where T<:Real
    if (interior_knots != nothing && boundary_knots != nothing)
        # pass
    elseif (knots != nothing)
        boundary_knots = extrema(knots)
        interior_knots = length(knots)==2 ? nothing : knots[2:(length(knots)-1)]
    else
        if (boundary_knots == nothing)
            boundary_knots = extrema(x)
        end
        iKnots = df - order + knots_offset + 1 - Int32(intercept)
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
`(x:: Array{T,1}; ders :: Int32 = 0)` for evaluation of `ders`
derivative for the splines at `x`.

The keyword arguments include one of:
1. `df`, possibly in combination with `intercept`
2. `boundary_knots` and `interior_knots`
3. `knots`

# Arguments
- `boundary_knots :: Union{Tuple{T,T},Nothing} = nothing`: boundary knots
- `interior_knots :: Union{Array{T,1},Nothing} = nothing`: interior knots
- `order :: Int32 = 4`: order of the spline
- `intercept :: Bool = false`: bool for whether to include an intercept
- `df :: Int32 = order - 1 + Int32(intercept)`: degrees of freedom
- `knots :: Union{Array{T,1}, Nothing} = nothing`: full set of knots
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
            order :: Int32 = 4,
            intercept :: Bool = false,
            df :: Int32 = order - 1 + Int32(intercept),
            knots :: Union{Array{T,1}, Nothing} = nothing,
            centre :: Union{T,Nothing} = nothing) where T<:Real
    (boundary_knots, interior_knots) =
        spline_args(x, boundary_knots=boundary_knots, interior_knots=interior_knots,
                    order=order, intercept=intercept,
                    df=df, knots=knots)
    spline = BSplineBasis(boundary_knots, interior_knots, order, intercept)
    function eval(x :: Array{T,1}; ders :: Int32 = 0)
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
- `order :: Int32 = 4`: order of the spline
- `intercept :: Bool = false`: bool for whether to include an intercept
- `df :: Int32 = order - 1 + Int32(intercept)`: degrees of freedom
- `knots :: Union{Array{T,1}, Nothing} = nothing`: full set of knots
- `centre :: Union{T,Nothing} = nothing)`: value to centre the splines
- `ders :: Int32 = 0`: derivatives of the splines

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
function bs(x :: Array{T,1}; ders :: Int32 = 0, kwargs...) where T<:Real
    bs_(x; kwargs...)(x, ders=ders)
end

"""
    ns_(x :: Array{T,1}; <keyword arguments>) where T<:Real

Calculate a basis for natural B-splines and return a function with signature
`(x:: Array{T,1}; ders :: Int32 = 0)` for evaluation of `ders`
derivative for the splines at `x`.

The keyword arguments include one of:
1. `df`, possibly in combination with `intercept`
2. `boundary_knots` and `interior_knots`
3. `knots`

# Arguments
- `boundary_knots :: Union{Tuple{T,T},Nothing} = nothing`: boundary knots
- `interior_knots :: Union{Array{T,1},Nothing} = nothing`: interior knots
- `order :: Int32 = 4`: order of the spline
- `intercept :: Bool = false`: bool for whether to include an intercept
- `df :: Int32 = order - 1 + Int32(intercept)`: degrees of freedom
- `knots :: Union{Array{T,1}, Nothing} = nothing`: full set of knots
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
            order :: Int32 = 4,
            intercept :: Bool = false,
            df :: Int32 = order - 3 + Int32(intercept),
            knots :: Union{Array{T,1}, Nothing} = nothing,
            centre :: Union{T,Nothing} = nothing) where T<:Real
    (boundary_knots, interior_knots) =
        spline_args(x, boundary_knots=boundary_knots, interior_knots=interior_knots,
                    order=order, intercept=intercept, df=df, knots=knots, knots_offset=2)
    spline = NSplineBasis(boundary_knots, interior_knots, order, intercept)
    function eval(x :: Array{T,1}; ders :: Int32 = 0)
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
- `order :: Int32 = 4`: order of the spline
- `intercept :: Bool = false`: bool for whether to include an intercept
- `df :: Int32 = order - 1 + Int32(intercept)`: degrees of freedom
- `knots :: Union{Array{T,1}, Nothing} = nothing`: full set of knots
- `centre :: Union{T,Nothing} = nothing)`: value to centre the splines
- `ders :: Int32 = 0`: derivatives of the splines

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
function ns(x :: Array{T,1}; ders :: Int32 = 0, kwargs...) where T<:Real
    ns_(x; kwargs...)(x, ders=ders)
end

"""
    is_(x :: Array{T,1}; <keyword arguments>) where T<:Real

Calculate a basis for I-splines and return a function with signature
`(x:: Array{T,1}; ders :: Int32 = 0)` for evaluation of `ders`
derivative for the splines at `x`.

The keyword arguments include one of:
1. `df`, possibly in combination with `intercept`
2. `boundary_knots` and `interior_knots`
3. `knots`

# Arguments
- `boundary_knots :: Union{Tuple{T,T},Nothing} = nothing`: boundary knots
- `interior_knots :: Union{Array{T,1},Nothing} = nothing`: interior knots
- `order :: Int32 = 4`: order of the spline
- `intercept :: Bool = false`: bool for whether to include an intercept
- `df :: Int32 = order - 1 + Int32(intercept)`: degrees of freedom
- `knots :: Union{Array{T,1}, Nothing} = nothing`: full set of knots

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
            order :: Int32 = 4,
            intercept :: Bool = false,
            df :: Int32 = order - 1 + Int32(intercept),
            knots :: Union{Array{T,1}, Nothing} = nothing) where T<:Real
    (boundary_knots, interior_knots) =
        spline_args(x, boundary_knots=boundary_knots, order=order+1,
                    intercept=intercept, df=df+1, knots=knots)
    spline = BSplineBasis(boundary_knots, interior_knots, order+1, false)
    knots = parent(spline.spline_basis.knots)
    findj(x) = interior_knots==nothing ? (order+1) : searchsortedlast(knots,x)
    function eval(x :: Array{T,1}; ders :: Int32 = 0)
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
                if j < js-(order+1)
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
- `order :: Int32 = 4`: order of the spline
- `intercept :: Bool = false`: bool for whether to include an intercept
- `df :: Int32 = order - 1 + Int32(intercept)`: degrees of freedom
- `knots :: Union{Array{T,1}, Nothing} = nothing`: full set of knots
- `ders :: Int32 = 0`: derivatives of the splines

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
function is(x :: Array{T,1}; ders :: Int32 = 0, kwargs...) where T<:Real
    is_(x; kwargs...)(x, ders=ders)
end

"""
    ms_(x :: Array{T,1}; <keyword arguments>) where T<:Real

Calculate a basis for M-splines and return a function with signature
`(x:: Array{T,1}; ders :: Int32 = 0)` for evaluation of `ders`
derivative for the splines at `x`.

The keyword arguments include one of:
1. `df`, possibly in combination with `intercept`
2. `boundary_knots` and `interior_knots`
3. `knots`

# Arguments
- `boundary_knots :: Union{Tuple{T,T},Nothing} = nothing`: boundary knots
- `interior_knots :: Union{Array{T,1},Nothing} = nothing`: interior knots
- `order :: Int32 = 4`: order of the spline
- `intercept :: Bool = false`: bool for whether to include an intercept
- `df :: Int32 = order - 1 + Int32(intercept)`: degrees of freedom
- `knots :: Union{Array{T,1}, Nothing} = nothing`: full set of knots
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
            order :: Int32 = 4,
            intercept :: Bool = false,
            df :: Int32 = order - 1 + Int32(intercept),
            knots :: Union{Array{T,1}, Nothing} = nothing,
            centre :: Union{T,Nothing} = nothing) where T<:Real
    (boundary_knots, interior_knots) =
        spline_args(x, boundary_knots=boundary_knots, interior_knots=interior_knots,
                    order=order, intercept=intercept, df=df, knots=knots)
    spline = BSplineBasis(boundary_knots, interior_knots, order, intercept)
    function eval(x :: Array{T,1}; ders :: Int32 = 0)
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
- `order :: Int32 = 4`: order of the spline
- `intercept :: Bool = false`: bool for whether to include an intercept
- `df :: Int32 = order - 1 + Int32(intercept)`: degrees of freedom
- `knots :: Union{Array{T,1}, Nothing} = nothing`: full set of knots
- `centre :: Union{T,Nothing} = nothing)`: value to centre the splines
- `ders :: Int32 = 0`: derivatives of the splines

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
function ms(x :: Array{T,1}; ders :: Int32 = 0, kwargs...) where T<:Real
    ms_(x; kwargs...)(x, ders=ders)
end

end # module
