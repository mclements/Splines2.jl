module Splines2StatsModelsExt

using Splines2
using StatsModels

using StatsModels

Splines2.ns(x, df::Integer) = Splines2.ns(x; df, intercept=false)

const SPLINE_CONTEXT = Any

mutable struct NSplineTerm{T,D} <: AbstractTerm
    const term::T
    const df::D
    basis::Union{NSplineBasis,Nothing}
end

Base.show(io::IO, p::NSplineTerm) = print(io, "ns($(p.term), $(p.df))")

function StatsModels.apply_schema(t::FunctionTerm{typeof(ns)},
                                  sch::StatsModels.Schema,
                                  Mod::Type{<:SPLINE_CONTEXT})
    length(t.args) == 2 || error("Incorrect number of arguments to ns")
    return apply_schema(NSplineTerm(t.args[1], t.args[2], nothing), sch, Mod)
end

function StatsModels.apply_schema(t::NSplineTerm,
                                  sch::StatsModels.Schema,
                                  Mod::Type{<:SPLINE_CONTEXT})
    term = apply_schema(t.term, sch, Mod)
    isa(term, ContinuousTerm) ||
        throw(ArgumentError("NSplineTerm only works with continuous terms (got $term)"))
    isa(t.df, ConstantTerm) ||
        throw(ArgumentError("NSplineTerm df must be a number (got $(t.df))"))
    return NSplineTerm(term, t.df.n, t.basis)
end

function StatsModels.modelcols(p::NSplineTerm{<:ContinuousTerm,<:Integer}, d::NamedTuple)
    col = modelcols(p.term, d)
    if isnothing(p.basis)
        p.basis = NSplineBasis(col; p.df)
    end
    return p.basis(col)
end

StatsModels.terms(p::NSplineTerm) = terms(p.term)
StatsModels.termvars(p::NSplineTerm) = StatsModels.termvars(p.term)
StatsModels.width(p::NSplineTerm) = 1
function StatsModels.coefnames(p::NSplineTerm)
    return string.("ns(", coefnames(p.term), ", ", 1:(p.df), ")")
end

end # module
