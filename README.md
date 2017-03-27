# A Module for Sparse Grid Discretization using Discontinuous Galerkin Bases

[![Build Status](https://travis-ci.org/ABAtanasov/GalerkinSparseGrids.jl.svg?branch=master)](https://travis-ci.org/ABAtanasov/GalerkinSparseGrids.jl)
[![codecov](https://codecov.io/gh/ABAtanasov/GalerkinSparseGrids.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ABAtanasov/GalerkinSparseGrids.jl)

This Julia Language Package is intended for accurately and efficiently solving hyperbolic partial differential equations in higher dimensions, where the curse of dimensionality restricts the computational feasibility of discretization of space using regular grids. Instead, we employ the sparse grid construction as in Bungartz & Griebel:
http://wissrech.ins.uni-bonn.de/research/pub/griebel/sparsegrids.pdf.

This technique in particular allows for an efficient numerical solution of Einstein's equations in full 3+1 dimensional space, as well as for tackling other systems in high-dimensional condensed matter calculations.

The sparse grid methods of this package can be extended beyond numerical relativity to many areas in high-dimensional data science and dynamics.

## Installing

Prerequisites for using this package are:

	Julia-0.5 (Latest Release)
	ODE.jl	(For running the timesteps to solve the sparse dynamical evolution)
	Cubature.jl (working to remove this prerequisite)

They can be pulled from <https://github.com/JuliaLang/ODE.jl> and <https://github.com/stevengj/Cubature.jl>, respectively.

Within Julia, use the package manager to write
`Pkg.add("GalerkinSparseGrids")` to locally install this package. 

The latest version is available to be pulled from <https://github.com/ABAtanasov/GalerkinSparseGrids.jl>. You can access it by running `git pull https://github.com/ABAtanasov/GalerkinSparseGrids.jl master` from the appropriate package directory.

## Functionality

This package allows for the efficient interpolation, differentiation, and time evolution of high-dimensional datasets.


### Summary of the DG Basis

In the 1-dimensional setting, the standard DG basis at order `k` and resolution `n` is equivalent to subdividing the axis in to 2^n sub-intervals and interpolating the function by polynomials of order strictly less than `k` on each of these sub-intervals. This gives a vector space of size `k 2^n`.

To make this compatible with sparse grids in higher dimensions, we find a "multi-resolution" basis for this vector space. That is, our basis functions `v(level, place, f_number)` are indexed by a level ranging from `0` to `n` and then a place ranging from `0` to `2^l - 1` so that `v` is supported on (p 2^(-l), (p+1) 2^(-l)). Lastly, `f_number` ranges over the `k` polynomials on each such interval. These `v` are made orthonormal to one another by a Grahm-Schmidt process detailed in Alpert and implemented in `DG_Basis.jl`

In D-dimensions, we use the standard tensor-products construction on these 1-D basis functions to obtain a basis of size `(k 2^n)^D`. For a default full grid, this vector space is equivalent to subdividing each axis of this domain into 2^n sub-intervals, for a total of `2^(n D)` hypercubes, and interpolating the function on each of these subdomains by multivariate polynomials of degree strictly less than k. The basis functions V(level, place, f_number) now have each of their arguments as a D-vector. 

Because this scheme becomes rapidly and prohibitively expensive in high dimensions, we apply the well-known sparse grid cutoff to exclude all basis functions with levels having 1-norm greater than n. This reduces us from O(k^D 2^(n D)) to O(k^D 2^n n^(D-1)).

### Interpolation:

To interpolate a given function f on the D-dimensional hypercube [0,1]^D, the `coeffs_DG(D::Int, k::Int, n::Int, f::Function; scheme="sparse")` method is used. 

For example, to interpolate a 2-D function at order k = 3 and n = 5, we would do

	D = 2; k = 3; n = 5;
	f = x->sin(2*pi*x[1])*sin(2*pi*x[2])
	full_coeffs   = coeffs_DG(D, k, n, f; scheme="full"  )
	sparse_coeffs = coeffs_DG(D, k, n, f; scheme="sparse")

This gives a dictionary indexed by multi-levels of type CartesianIndex{2}. For a given multi-level, any dictionary entry is an array of arrays indexed by the multi-place and multi-fnumber that holds the coefficient data.

In either case, the interpolation can be evaluated at a given point xs (given by an array of length D) by 

	full_val   = reconstruct_DG(full_coeffs  , xs)
	sparse_val = reconstruct_DG(sparse_coeffs, xs)
	
No scheme, D, k, or n, need be specified, as all of these can be deduced from the dictionary itself. 

The default scheme is always "sparse".

### Solving Hyperbolic PDEs

All solvers can be found in the `PDEs.jl` script. The simplest one solves the wave equation in 1-D using the DG position basis 
	
	wave_evolve_1D(k::Int, n::Int, f0::Function, v0::Function, time0::Real, time1::Real; base = "hier", order = "45"). 
	
Here, k and n are as before, f0 is the initial condition for position, v0 is the initial condition for velocity, time0 and time1 are the initial and final times of evolution. For 1D, we can also use the "pos" position basis (rather than the multi-resolution), and order "78" ode solvers. In either case, `ode45/ode78` from `ODE.jl` are respectively used to solve the differential equation.
    
For example, for a standing wave solution from t_0 to t_1 seconds with k polynomials on each division, up to hierarchical level max_level, using ode78:

	f0 = x->sin(2*pi*x[1])
	v0 = x->0
    hier_soln = wave_evolve_1D(k, max_level, f0, v0, t_0, t_1; order="78")

For higher dimensions, we use 

	wave_evolve(D::Int, k::Int, n::Int, f0::Function, v0::Function, time0::Real, time1::Real; order = "45", scheme="sparse").

For example, to solve a standing wave in 2-D from t_0 to t_1 seconds at k polynomials on each division, sparse depth n and using ode78:

	f0 = x->sin(2*pi*x[1])*sin(2*pi*x[2])
	v0 = 0
	sparse_soln = wave_evolve(2, k, n, f0, v0, t_0, t_1; order="78")

At n = 1, this is the same as `wave_evolve_1D` in the "hier" basis.  

Although the option scheme="full" is implemented for higher dimensions, it is almost never computationally feasible, even in dimensions as low as D=3.

### Coefficient Vectors

It is best to work directly with vectors of coefficients rather than dictionary or tree-like structures when defining linear operators (like differentiation) acting on the interpolations. This way, derivatives are represented by sparse matrices that can act on the coefficient vectors with very fast performance. This is ideal for time evolution.

For this reason we have the `vcoeffs_DG` method that works exactly the same as their non-vector relatives, but return a coefficient vector instead of a dictionary. 

We can convert between coefficient vectors and corresponding dictionaries as follows (in 2D with k = 3, n = 5):

	D = 2; k = 3; n = 5;
	f = x->sin(4*x[1]+x[2])
	
	full_vcoeffs = vcoeffs_DG(D, k, n, f; scheme="full")
	full_dict = V2D(D, k, n, full_vcoeffs; scheme="full")
	# Note this also means full_vcoeffs = D2V(D, k, n, full_dict; scheme="full") 

	sparse_vcoeffs = vcoeffs_DG(D, k, n, f; scheme="sparse")
	sparse_dict = V2D(3,sparse_vcoeffs,5,2)
	# Note this also means sparse_vcoeffs = D2V(D, k, n, sparse_dict; scheme="sparse") 

In addition, we can generate a lookup by using the `D2Vref` and `V2Dref`that takes us from a dictionary-like level-place-fnumber scheme to a specific index in a coefficient vector and vice versa. This is useful when trying to understand specific values in a vector of coefficients. 

If we have the vector `sparse_vcoeffs` as above and wanted to understand what (level, place, fnumber) corresponds to the value of sparse_vcoeffs[10], the code would be

	VD = V2Dref(D, k, n; scheme="sparse")
	VD[10]

Conversely, If we wanted to see what index we should look at to get the coefficient for (level, place, fnumber) = ((4,4), (1,2), (2, 3)), the code would be

	level 	= CartesianIndex{2}((4, 4))
	place 	= CartesianIndex{2}((1, 2))
	fnumber = CartesianIndex{2}((2, 3))
	DV = D2Vref(D, k, n; scheme="sparse")
	DV[(level, place, fnumber)]

### Differentiation

The higher dimensional full and sparse derivative operators are implemented in `Multidim_Derivative.jl`. Here, `D`, `k`, and `n` are as before, and `i` specifies the axis along which the derivative is taken.

	D_matrix(D::Int, i::Int, k::Int, n::Int; scheme="sparse")
	
We also have the gradient vector that is `D_matrix` over all i, and the Laplacian:

	grad_matrix(D::Int, k::Int, n::Int; scheme="sparse")
	laplacian_matrix(D::Int, k::Int, n::Int; scheme="sparse")

All these matrices act directly on the coefficient vectors obtained from `vcoeffs_DG`. For now, we assume periodic boundary conditions. It is not too difficult to generalize away from a periodic boundary. 

In the discontinuous Galerkin method, it is customary to formulate the derivative in the weak sense. This is how the method is implemented. As is standard, this is done through an integration by parts, leading to the derivative matrix having two summands: one being a boundary term and the other a derivative term on the interior. 

## Other Bases

### Canonical "Hat" Basis

To interpolate a multivariate function in a standard "position" basis, use the `standard_coefficients(f::Function, ls::NTuple{D,Int})`. This takes a scalar function that takes an `D`-dimensional array input `xs` and interpolates it to resolution `2^l[i]` along the `i`th coordinate direction.

From this interpolation, we can reconstruct the function at a point `xs` from the basis coefficients by using the `standard_reconstruct(coefficients::AbstractArray, ls::NTuple{D,Int}, xs::NTuple{D,T})` function.

To interpolate a multivariate function in the multiresolution (aka hierarchical) scheme (c.f. Gerstner & Griebel below) use `coeffs_hat(D::Int, n::Int, f::Function)`. This gives the exact same interpolation but within the multiresolution basis. The coefficients are given as a dictionary of CartesianIndex types. The keys of the dictionary correspond to multi-levels and corresponding CartesianIndex entries correspond to the set of places at a specific multilevel.

Similarly, we can interpolate in a sparse basis (equivalent to the multiresolution basis in 1D) using `sparse_coefficients(D::Int, n::Int, f::Function)` where `n` is the max 1-norm of the multi-level and `D` is the dimension of the space. For example in 2-D:

	scoeffs = coeffs_hat(2, 5, 2x->sin(4*x[1]+x[2]))

For either the hierarchical and sparse basis coefficients, we can reconstruct the interpolation at a point xs using the `reconstruct_hat(coefficients::Dict{CartesianIndex{D},Array{Float64,D}}, x::NTuple{D,T})` function. For example we could take the prior result and evaluate at the point (.4,.6) by

	reconstruct_hat(scoeffs, (.4,.6))

## Construction of the DG Basis

The DG position basis in consists of `k` Legendre polynomials (from degree 0 to `k-1`) supported on subintervals of the form `(i 2^-level, (i+1) 2^-level)`. By construction, this gives an orthonormal basis. For this standard 1-D DG basis, the method `get_vcoeffs(k::Int, level::Int, f::Function)` gives an appropriate array of coefficients that can be reconstructed at a point `x` using `reconstruct_vcoeffs(k,level,coefficients,x)`. 

The corresponding multiresolution basis formed by a series of orthogonalizations implemented in `DG_Functions.jl`. This basis spans the same space as the position basis, but makes use of a discontinuous basis (c.f. `1D_DG_Functions.jl`) to achieve the hierarchical structure.

(More to come on orthogonalization schemes)


## Future updates

The family of PDEs that can be solved is currently being extended to include a wide variety of "Einstein-like" equations.

Future updates include adaptivity as well as parallelizability to these methods. 

## References

A Sparse Grid Discontinuous Galerkin Method for high-dimensional transport equations
http://arxiv.org/abs/1602.02124

(also c.f. for elliptic equations: http://arxiv.org/abs/1508.07781)

Gerstner & Griebel:
A survey of sparse grid discretization and applciations to quadrature, data compression, and dynamical systems: 
http://wissrech.ins.uni-bonn.de/research/pub/gerstner/sparsegrids.pdf

Good introduction to the Interior Penality Discontinuous Galerkin (IPDG):
http://www.cs.elte.hu/~izsakf/otka/presentation_dg_23dim_en.pdf
