module Splines2

export ns, NSplineBasis, bs, BSplineBasis, is, ISplineBasis, ms, MSplineBasis, basis

using LinearAlgebra
using Statistics

abstract type AbstractSplineBasis{T<:Real} end

include("basis.jl")
include("bspline.jl")
include("nspline.jl")
include("ispline.jl")
include("mspline.jl")
include("utils.jl")

end # module
