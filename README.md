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

Similarly, we can interpolate in a sparse basis (equivalent to the multiresolution basis in 1D) using `sparse_coefficients(f::Function, n::Int, D::Int)` where `n` is the max 1-norm of the multi-level and `D` is the dimension of the space. For example in 2-D:

    scoeffs = sparse_coefficients(x->sin(4*x[1]+x[2]),5,2)


For either the hierarchical and sparse basis coefficients, we can reconstruct the interpolation at a point xs using the `reconstruct(coefficients::Dict{CartesianIndex{D},Array{Float64,D}}, x::NTuple{D,T})` function. For example we could take the prior result and evaluate at the point (.4,.6) by

    reconstruct(scoeffs, (.4,.6))

### DG Basis

The DG position basis in consists of `k` Legendre polynomials (from degree 0 to `k-1`) supported on subintervals of the form `(i 2^-level, (i+1) 2^-level)`. By construction, this gives an orthonormal basis. For this standard 1-D DG basis, the method `get_vcoeffs(k::Int, level::Int, f::Function)` gives an appropriate array of coefficients that can be reconstructed at a point `x` using `reconstruct_vcoeffs(k,level,coefficients,x)`. 

The corresponding multiresolution basis formed by a series of orthogonalizations implemented in `DG_Functions.jl`. This basis spans the same space as the position basis, but makes use of discontinuous basis (c.f. `Specific_DG_Functions.jl`) to achieve the hierarchical structure.

We can interpolate multi-dimensional functions using the DG basis by using the methods `hier_coefficients_DG(k::Int, f::Function, ls::NTuple{D,Int})` and `hier_coefficients_DG(k::Int, f::Function, n::Int, D::Int)` with the same syntax as the corresponding non-DG methods. Similarly, to reconstruct a function at a point `xs` from a coefficient dictionary we use `reconstruct_DG(k,coefficients,xs)`.

### Coefficient Vectors

It is best to work directly with vectors of coefficients rather than dictionary or tree-like structures when defining differentiation operators. This way, derivatives are represented by sparse matrices that can act on the coefficient vectors with very fast performance. This is ideal for time evolution.

For this reason we have the `vhier_coefficients_DG` and `vsparse_coefficients_DG` methods that work exactly the same as their non-vector relatives, but return a coefficient vector instead of a dictionary. 

We can convert between vectors and appropriate dictionaries as follows (for k = 3, l = (5,5) or n=5):

    full_vcoeffs = vhier_coefficients_DG(3,x->sin(4*x[1]+x[2]),(5,5))
    full_dict = Full_V2D(3,full_vcoeffs,(5,5))
    # Note this also means full_vcoeffs = Full_D2V(3,full_dict,(5,5)) 

    sparse_vcoeffs = vsparse_coefficients_DG(3,x->sin(4*x[1]+x[2]),5,2)
    sparse_dict = Sparse_V2D(3,sparse_vcoeffs,5,2)
    # Note this also means sparse_vcoeffs = Sparse_D2V(3,sparse_dict,5) 

In addition, we can generate a lookup by using the `full_referenceD2V`/`V2D` and `sparse_referenceD2V`/`V2D` that takes us from a dictionary-like level-place-fnumber scheme to a specific place in a coefficient vector and vice versa. 

### Differentiation Schemes

In the discontinuous Galerkin method, it is customary to formulate the derivative in the weak sense. Here, we do the same. As is standard, we perform an integration by parts, leading to the derivative matrix having two summands: one a boundary term and the other a derivative term on the interior. It is the boundary term that we must be careful with. Although using a Lax-Friedrichs flux is a common choice, we avoid this (effectively choosing alpha = 0). 

For now, we assume periodic boundary conditions. It is not too difficult to generalize away from a periodic boundary. 

The derivative matrix for the position basis is `periodic_pos_DLF_Matrix(alpha::Real, k::Int, max_level::Int)` while the matrix for the hierarchical basis is obtained from this by conjugation and can be called by `periodic_hier_DLF_Matrix(alpha::Real, k::Int, max_level::Int)`. Here DLF stands for Derivative with Lax-Friedrichs flux. 

It is advised, for now, to have alpha = 0 (so no Lax-Friedrichs flux is involved).

### Solving Hyperbolic PDEs

All solvers can be found in the `PDEs.jl` script. The simplest one solves the wave equation in 1-D using the DG position basis `pos_wave_equation45(f0::Function, v0::Function, k::Int,level::Int, time0::Real, time1::Real)`. The returned value is the same as that returned by `ode45` in `ODE.jl`
    
For example, a standing wave solution from 0 to 1 seconds with k=3, l=5 is given by:

    pos_soln = pos_wave_equation45(x->sin(2*pi*x[1]), x->0, 3, 5, 0, 1)

The same syntax applies to solving the problem in the hierarchical 1-D basis using `hier_wave_equation45`. This is equivalent to the position basis solution. 

For higher dimensions, we use sparse grids to yield the methods `sparse_wave_equation45` and `sparse_wave_equation78`. For example, to solve a standing wave in 2-D from 0 to 1 seconds at k=4, n = 6 we can use

    sparse_soln = sparse_wave_equation45(x->sin(2*pi*x[1])*sin(2*pi*x[2]), x-> 0, 4, 6, 2, 0, 1)

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
