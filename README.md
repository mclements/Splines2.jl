# Splines2.jl package for regression splines

A [Julia](https://julialang.org/) package for regression splines. The package currently includes B-splines, natural B-splines, M-splines and I-splines.

## News
### Version 0.2.0:
- Mainly bug fixes.
- A change of behaviour for `Splines2.is_` and `Splines2.is`: `intercept=true` will include a columns of ones, while the default `intercept=false` will keep _all_ of the spline terms, but exclude the column of ones. This behaviour is different to the `splines2` package in R, which will give all of the spline terms for `intercept=TRUE` and drop the first spline term for `intercept=FALSE`.

## Installation

The package is registered on JuliaHub. For installation:

``` julia
using Pkg; Pkg.add("Splines2")
```

## Usage

Exported functions include `Splines2.bs`, `Splines2.ns`, `Splines2.ms` and `Splines2.is`, which provide evaluating spline bases for B-splines, natural B-splines, M-splines and I-splines, respectively. These functions take an `::Array{<:Real,1}` argument and some design information and return the given spline basis. 

### Documentation for `Splines2.bs`

``` julia
bs(x :: Array{T,1}; <keyword arguments>) where T<:Real
```

Calculate a basis for B-splines. 

The keyword arguments include one of:
1. `df`, possibly in combination with `intercept`
2. `boundary_knots` and `interior_knots`
3. `knots`

#### Arguments
- `boundary_knots :: Union{Tuple{T,T},Nothing} = nothing`: boundary knots
- `interior_knots :: Union{Array{T,1},Nothing} = nothing`: interior knots
- `order :: Int32 = 4`: order of the spline
- `intercept :: Bool = false`: bool for whether to include an intercept
- `df :: Int32 = order - 1 + Int32(intercept)`: degrees of freedom
- `knots :: Union{Array{T,1}, Nothing} = nothing`: full set of knots (excluding repeats)
- `centre :: Union{T,Nothing} = nothing)`: value to centre the splines
- `ders :: Int32 = 0`: derivatives of the splines



### Documentation for `Splines2.bs_`

``` julia
bs_(x :: Array{T,1}; <keyword arguments>) where T<:Real
```

Calculate a basis for B-splines and return a function with signature
`(x:: Array{T,1}; ders :: Int32 = 0)` for evaluation of `ders`
derivative for the splines at `x`.

The keyword arguments include one of:
1. `df`, possibly in combination with `intercept`
2. `boundary_knots` and `interior_knots`
3. `knots`

#### Arguments
- `boundary_knots :: Union{Tuple{T,T},Nothing} = nothing`: boundary knots
- `interior_knots :: Union{Array{T,1},Nothing} = nothing`: interior knots
- `order :: Int32 = 4`: order of the spline
- `intercept :: Bool = false`: bool for whether to include an intercept
- `df :: Int32 = order - 1 + Int32(intercept)`: degrees of freedom
- `knots :: Union{Array{T,1}, Nothing} = nothing`: full set of knots (excluding repeats)
- `centre :: Union{T,Nothing} = nothing)`: value to centre the splines

The documentation for the other bases are similar, except that the I-splines do not include the `centre` argument.

## Examples

Some short examples are given below.

``` julia
julia> using Splines2
julia> x = collect(0.0:0.1:1.0);
julia> bs(x, df=3)

11×3 Array{Float64,2}:
 0.0    0.0    0.0  
 0.243  0.027  0.001
 0.384  0.096  0.008
 0.441  0.189  0.027
 0.432  0.288  0.064
 0.375  0.375  0.125
 0.288  0.432  0.216
 0.189  0.441  0.343
 0.096  0.384  0.512
 0.027  0.243  0.729
 0.0    0.0    1.0
 
julia> ns(x, boundary_knots=(0.0,1.0), interior_knots=[0.2])
	
11×2 Array{Float64,2}:
 0.0        0.0      
 0.196457  -0.106365 
 0.363908  -0.179949 
 0.479393  -0.194802 
 0.544119  -0.152288 
 0.565337  -0.0606039
 0.550299   0.072056 
 0.506256   0.237496 
 0.44046    0.427522 
 0.360161   0.633938 
 0.272611   0.84855
 
 julia> ms(x, knots=[0.0,0.4,1.0], centre=0.4)

11×4 Array{Float64,2}:
 -1.44      -1.92       -0.64      0.0      
  0.6075    -1.665      -0.63      0.0      
  1.14      -1.08       -0.56      0.0      
  0.7425    -0.435      -0.37      0.0      
  0.0        0.0         0.0       0.0      
 -0.606667   0.0244444   0.563704  0.0308642
 -1.01333   -0.284444    1.14963   0.246914 
 -1.26      -0.78        1.54      0.833333 
 -1.38667   -1.31556     1.51704   1.97531  
 -1.43333   -1.74444     0.862963  3.85802  
 -1.44      -1.92       -0.64      6.66667
 
 julia> is(x, df=4)

11×4 Array{Float64,2}:
 0.0     0.0     0.0     0.0   
 0.3439  0.0523  0.0037  0.0001
 0.5904  0.1808  0.0272  0.0016
 0.7599  0.3483  0.0837  0.0081
 0.8704  0.5248  0.1792  0.0256
 0.9375  0.6875  0.3125  0.0625
 0.9744  0.8208  0.4752  0.1296
 0.9919  0.9163  0.6517  0.2401
 0.9984  0.9728  0.8192  0.4096
 0.9999  0.9963  0.9477  0.6561
 1.0     1.0     1.0     1.0   
```

We also provide functions that return a function for evaluating spline bases with a function signature `(x::Array{T<:Real,1}; ders::Int32 = 0)`. These are useful for "safe" predictions in regression modelling. As an example:

``` julia
julia> using Splines2, GLM, Random
julia> Random.seed!(12345);
julia> x = collect(range(0.0, length=301, stop=2.0*pi));
julia> y = sin.(x)+randn(length(x)); 
julia> ns1 = Splines2.ns_(x,df=5,intercept=true); # this is a function
julia> X = ns1(x);
julia> fit1 = lm(X,y)

LinearModel{GLM.LmResp{Array{Float64,1}},GLM.DensePredChol{Float64,LinearAlgebra.Cholesky{Float64,Array{Float64,2}}}}:

Coefficients:
────────────────────────────────────────────────────────────────────
     Estimate  Std. Error    t value  Pr(>|t|)  Lower 95%  Upper 95%
────────────────────────────────────────────────────────────────────
x1   1.23751     0.269035   4.59981     <1e-5    0.708047   1.76698 
x2   0.12448     0.249256   0.499407    0.6179  -0.366058   0.615018
x3  -1.89278     0.256808  -7.37043     <1e-11  -2.39819   -1.38738 
x4   0.187169    0.22469    0.833012    0.4055  -0.255023   0.629361
x5  -0.240554    0.254986  -0.943404    0.3462  -0.742369   0.26126 
────────────────────────────────────────────────────────────────────

julia> newx = collect(0.0:0.5:3.5);
julia> predict(fit1, ns1(newx)) # safe predictions

8-element Array{Float64,1}:
  0.2982757838333453 
  0.6021897830602807 
  0.8365641389496451 
  0.9318592081638681 
  0.8310124040845238 
  0.5536590079608558 
  0.14855743047881534
 -0.3299373638222967 
```


## Using `Splines2` with `@formula`

We provide code below for using the `Splines2` package with `@formula`. Note that these do *not* provide "safe" predictions.

``` julia
using StatsModels
ns(x,df) = Splines2.ns(x,df=df,intercept=true) # assumes intercept
const NSPLINE_CONTEXT = Any
struct NSplineTerm{T,D} <: AbstractTerm
    term::T
    df::D
end
Base.show(io::IO, p::NSplineTerm) = print(io, "ns($(p.term), $(p.df))")
function StatsModels.apply_schema(t::FunctionTerm{typeof(ns)},
                                  sch::StatsModels.Schema,
                                  Mod::Type{<:NSPLINE_CONTEXT})
    apply_schema(NSplineTerm(t.args_parsed...), sch, Mod)
end
function StatsModels.apply_schema(t::NSplineTerm,
                                  sch::StatsModels.Schema,
                                  Mod::Type{<:NSPLINE_CONTEXT})
    term = apply_schema(t.term, sch, Mod)
    isa(term, ContinuousTerm) ||
        throw(ArgumentError("NSplineTerm only works with continuous terms (got $term)"))
    isa(t.df, ConstantTerm) ||
        throw(ArgumentError("NSplineTerm df must be a number (got $(t.df))"))
    NSplineTerm(term, t.df.n)
end
function StatsModels.modelcols(p::NSplineTerm, d::NamedTuple)
    col = modelcols(p.term, d)
    Splines2.ns(col, df=p.df,intercept=true)
end
StatsModels.terms(p::NSplineTerm) = terms(p.term)
StatsModels.termvars(p::NSplineTerm) = StatsModels.termvars(p.term)
StatsModels.width(p::NSplineTerm) = 1
StatsModels.coefnames(p::NSplineTerm) = "ns(" .* coefnames(p.term) .* "," .* string.(1:p.df) .* ")"
```

To show that this is not safe:

``` julia
julia> using DataFrames
julia> d = DataFrames.DataFrame(x=x,y=y);
julia> fit2 = lm(@formula(y~ns(x,5)+0),d) # equivalent to fit1 with nicer labels

StatsModels.TableRegressionModel{LinearModel{GLM.LmResp{Array{Float64,1}},GLM.DensePredChol{Float64,LinearAlgebra.Cholesky{Float64,Array{Float64,2}}}},Array{Float64,2}}

y ~ 0 + ns(x, 5)

Coefficients:
─────────────────────────────────────────────────────────────────────────
          Estimate  Std. Error    t value  Pr(>|t|)  Lower 95%  Upper 95%
─────────────────────────────────────────────────────────────────────────
ns(x,1)   1.23751     0.269035   4.59981     <1e-5    0.708047   1.76698 
ns(x,2)   0.12448     0.249256   0.499407    0.6179  -0.366058   0.615018
ns(x,3)  -1.89278     0.256808  -7.37043     <1e-11  -2.39819   -1.38738 
ns(x,4)   0.187169    0.22469    0.833012    0.4055  -0.255023   0.629361
ns(x,5)  -0.240554    0.254986  -0.943404    0.3462  -0.742369   0.26126 
─────────────────────────────────────────────────────────────────────────

julia> predict(fit2, DataFrames.DataFrame(x=newx)) # unsafe predictions!

8-element Array{Union{Missing, Float64},1}:
  0.29827578383334535
  0.7976143687096604 
  0.8964195213823501 
  0.40991870738161984
 -0.4167421148184624 
 -1.0400611367418444 
 -0.7710405835443831 
  0.1787886772299305
```

For further details, see the discussion [here](https://discourse.julialang.org/t/safe-predictions-using-formula-and-regression-splines/33057).
