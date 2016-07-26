# A Module for Sparse Grid Discretization using Discontinuous Galerkin Bases

[![Build Status](https://travis-ci.org/AlexAtanasov14/GalerkinSparseGrids.jl.svg?branch=master)](https://travis-ci.org/AlexAtanasov14/GalerkinSparseGrids.jl)
[![codecov](https://codecov.io/gh/AlexAtanasov14/GalerkinSparseGrids.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/AlexAtanasov14/GalerkinSparseGrids.jl)

This Julia Language Package is intended for accurately and efficiently sovling hyperbolic partial differential equations in higher dimensions, where the curse of dimensionality restricts the computational feasibility of discretization using regular grids. Instead, we employ the sparse grid construction as in Bungartz & Griebel:
http://wissrech.ins.uni-bonn.de/research/pub/griebel/sparsegrids.pdf

This technique in particular allows for an efficient numerical solution of Einstein's equations in full 3+1 dimensional space, as well as for tackling other systems in high-dimensional condensed matter calculations. 

## Installing

Within Julia, use the package manager to write
`Pkg.add("GalerkinSparseGrids")` to locally install this package. 

The latest version is available from <https://github.com/AlexAtanasov14/GalerkinSparseGrids.jl>. You can access it by running `git pull https://github.com/AlexAtanasov14/GalerkinSparseGrids.jl master` from the appropriate package directory.

## Functionality

This package allows for the quadrature, interpolation, and time evolution of high-dimensional datasets. 


## Future updates

We intend to include adaptivity as well as parallelizability to these methods. 

## References

A sparse grid discontinuous Galerkin method for high-dimensional transport equations
http://arxiv.org/abs/1602.02124

(also c.f. for elliptic equations: http://arxiv.org/abs/1508.07781)

A survey of sparse grid discretization and applciations to quadrature, data compression, and dynamical systems: 
http://wissrech.ins.uni-bonn.de/research/pub/gerstner/sparsegrids.pdf

Good introduction to the Interior Penality Discontinuous Galerkin (IPDG):
http://www.cs.elte.hu/~izsakf/otka/presentation_dg_23dim_en.pdf
