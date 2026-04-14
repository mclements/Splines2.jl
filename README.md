# Splines2.jl package for regression splines

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Stable Docs][docs-stable-img]][docs-stable-url]
[![Dev Docs][docs-dev-img]][docs-dev-url]
[![codecov](https://codecov.io/gh/mclements/Splines2.jl/graph/badge.svg?token=4559XBNT50)](https://codecov.io/gh/mclements/Splines2.jl)
[![Code Style: YAS](https://img.shields.io/badge/code%20style-yas-1fdcb2.svg)](https://github.com/jrevels/YASGuide)

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://mclements.github.io/Splines2.jl/dev

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://mclements.github.io/Splines2.jl/stable

A [Julia](https://julialang.org/) package for regression splines. The package currently includes B-splines, natural B-splines, M-splines and I-splines.

## News

### Version 0.3.0
- Support for safe predictions in `StatsModels.@formula` for natural splines.
- Substantial refactoring of internal code to use callable types.

### Version 0.2.0:
- Mainly bug fixes.
- A change of behaviour for `Splines2.is_` and `Splines2.is`: `intercept=true` will include a columns of ones, while the default `intercept=false` will keep _all_ of the spline terms, but exclude the column of ones. This behaviour is different to the `splines2` package in R, which will give all of the spline terms for `intercept=TRUE` and drop the first spline term for `intercept=FALSE`.

## Installation

The package is registered inn the Julia General registry. For installation:

``` julia
using Pkg; Pkg.add("Splines2")
```
