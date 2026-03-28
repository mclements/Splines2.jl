using Documenter
using Splines2

makedocs(;
         repo=Remotes.GitHub("palday", "Splines2.jl"),
         sitename="Splines2",
         doctest=true,
         checkdocs=:exports,
         warnonly=[:cross_references],
         format=Documenter.HTML(; edit_link="main"),
         pages=["index.md",
                "api.md"])

deploydocs(; repo="github.com/mclements/Splines2.jl.git",
           devbranch="main",
           push_preview=true)
