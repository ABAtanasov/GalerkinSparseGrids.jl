# A Module for Sparse Grid Discretization using Discontinuous Galerkin Bases

[![Build Status](https://travis-ci.org/ABAtanasov/GalerkinSparseGrids.jl.svg?branch=master)](https://travis-ci.org/ABAtanasov/GalerkinSparseGrids.jl)
[![codecov](https://codecov.io/gh/ABAtanasov/GalerkinSparseGrids.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ABAtanasov/GalerkinSparseGrids.jl)

This Julia Language Package is intended for accurately and efficiently sovling hyperbolic partial differential equations in higher dimensions, where the curse of dimensionality restricts the computational feasibility of discretization using regular grids. Instead, we employ the sparse grid construction as in Bungartz & Griebel:
http://wissrech.ins.uni-bonn.de/research/pub/griebel/sparsegrids.pdf

This technique in particular allows for an efficient numerical solution of Einstein's equations in full 3+1 dimensional space, as well as for tackling other systems in high-dimensional condensed matter calculations. 

## Installing

Within Julia, use the package manager to write
`Pkg.add("GalerkinSparseGrids")` to locally install this package. 

The latest version is available from <https://github.com/AlexAtanasov14/GalerkinSparseGrids.jl>. You can access it by running `git pull https://github.com/AlexAtanasov14/GalerkinSparseGrids.jl master` from the appropriate package directory.

## Functionality

This package allows for the quadrature, interpolation, and time evolution of high-dimensional datasets. 

### Canonical "Hat" Basis

To interpolate a multivariate function in a standard "position" basis, use the `standard_coefficients(f::Function, ls::NTuple{D,Int})`. This takes a scalar function that takes an `D`-dimensional array input `xs` and interpolates it to resolution `2^l[i]` along the `i`th coordinate direction.

From this interpolation, we can reconstruct the function at a point `xs` from the basis coefficients by using the `standard_reconstruct(coefficients::AbstractArray, ls::NTuple{D,Int}, xs::NTuple{D,T})` function.

To interpolate a multivariate function in the multiresolution (aka hierarchical) scheme (c.f. Gestner & Griebel below) use `hier_coefficients(f::Function, ls::NTuple{D,Int})`. This gives the exact same interpolation but within the multiresolution basis. The coefficients are given as a dictionary of CartesianIndex types. The keys of the dictionary correspond to multi-levels and corresponding CartesianIndex entries correspond to the set of places at a specific multilevel.

Similarly, we can interpolate in a sparse basis (equivalent to the multiresolution basis in 1D) using `hier_coefficients(f::Function, n::Int, D::Int)` where `n` is the max 1-norm of the multi-level and `D` is the dimension of the space. For example in 2-D:

    sparse_coeffs = sparse_coefficients(x->sin(4*x[1]+x[2]),l,2)


For either the hierarchical and sparse basis coefficients, we can reconstruct the interpolation at a point xs using the `reconstruct(coefficients::Dict{CartesianIndex{D},Array{Float64,D}}, x::NTuple{D,T})` function. For example we could take the prior result and evaluate at the point [.4,.6] by

    reconstruct(sparse_coeffs, [.4,.6])


### DG Basis

The DG piecewise polynomial basis is implemented in DG_Methods.jl


### Differentiation Schemes

### Solving Hyperbolic PDEs

## Future updates

The family of PDEs that can be solved is currently being extended to include a wide variety of "Einstein-like" equations.

We intend to include adaptivity as well as parallelizability to these methods. 

## References

A Sparse Grid Discontinuous Galerkin Method for high-dimensional transport equations
http://arxiv.org/abs/1602.02124

(also c.f. for elliptic equations: http://arxiv.org/abs/1508.07781)

Gerstner & Griebel:
A survey of sparse grid discretization and applciations to quadrature, data compression, and dynamical systems: 
http://wissrech.ins.uni-bonn.de/research/pub/gerstner/sparsegrids.pdf

Good introduction to the Interior Penality Discontinuous Galerkin (IPDG):
http://www.cs.elte.hu/~izsakf/otka/presentation_dg_23dim_en.pdf
