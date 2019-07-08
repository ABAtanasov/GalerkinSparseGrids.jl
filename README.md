# A Module for Sparse Grid Discretization using Discontinuous Galerkin Bases
[![Build Status](https://travis-ci.org/ABAtanasov/GalerkinSparseGrids.jl.svg?branch=master)](https://travis-ci.org/ABAtanasov/GalerkinSparseGrids.jl)
[![codecov](https://codecov.io/gh/ABAtanasov/GalerkinSparseGrids.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ABAtanasov/GalerkinSparseGrids.jl)


This [Julia language](https://julialang.org/) package is intended for accurately and efficiently solving hyperbolic partial differential equations in higher dimensions, where the curse of dimensionality restricts the computational feasibility of discretization of space using regular grid methods. Instead, we employ the sparse grid construction as in [Bungartz & Griebel](http://wissrech.ins.uni-bonn.de/research/pub/griebel/sparsegrids.pdf).

The ambitious long-term goal of this package is the efficient and accurate numerical solution of Einstein's equations in full 3+1 dimensional space under conditions of very low to no symmetry. Such simulations would be a great step forward in many diverse areas of strong gravity, including the modeling of dark matter, dark energy, and gravitational wave cosmology. We are not there yet, but we have already made progress in the simpler cases of transport and wave equation evolution and high-dimensional interpolation of initial conditions.

In order to make full anisotropic gravitational evolution a reality, several necessary additions are in the progress of being implemented:

* Handling a wider class of boundary value problems, accounting for decay conditions at infinity
* Handling one or more singularity points (e.g. for interacting black holes)
* Easy and general way to handle coordinate patching of the spacetime

For more information, see [Sparse Grid Discretizations based on a Discontinuous Galerkin Method](https://arxiv.org/abs/1710.09356)

The sparse grid methods of this package apply well beyond numerical relativity to many areas in high-dimensional dynamics and data science. The user is invited to experiment.

Authors: Alex Atanasov (atanasov@g.harvard.edu) & Erik Schnetter (eschnetter@perimeterinstitute.ca)

## Installing

Prerequisites for using this package are:

	Julia-1.0 (Latest Release)
	ODE.jl	(For running the timesteps to solve the sparse dynamical evolution)
	Cubature.jl (working to remove this prerequisite)
	StaticArrays (new, trying this out - could produce large speedup)

They can be added using the Julia package manager. Source documentation is [here](<https://github.com/JuliaLang/ODE.jl>) for ODE.jl and [here](<https://github.com/stevengj/Cubature.jl>) for Cubature.jl.

Within Julia, use the package manager to write `Pkg.clone("git://github.com/ABAtanasov/GalerkinSparseGrids.jl")` to locally install this package.

The latest version is available to be pulled from <http://github.com/ABAtanasov/GalerkinSparseGrids.jl>. You can access it by running `git pull https://github.com/ABAtanasov/GalerkinSparseGrids.jl master` from the appropriate package directory.

## Functionality

This package allows for the efficient interpolation, differentiation, and time evolution of high-dimensional datasets.


### Summary of the DG Basis

In the 1-dimensional setting, the standard DG basis at order `k` and resolution `n` is equivalent to subdividing the axis in to `2^n` sub-intervals and interpolating the function by polynomials of order strictly less than `k` on each of these sub-intervals. This gives a vector space of size `k 2^n`.

To make this compatible with sparse grids in higher dimensions, we find a "multi-resolution" basis for this vector space. That is, our basis functions `v(level, cell, mode)` are indexed by a level ranging from `0` to `n` and then a cell ranging from `0` to `2^l - 1` so that `v` is supported on (p 2^(-l), (p+1) 2^(-l)). Lastly, `mode` ranges over the `k` polynomials on each such interval. These `v` are made orthonormal to one another by a Grahm-Schmidt process detailed in Alpert and implemented in `DG_Basis.jl`

In D-dimensions, we use the standard tensor-products construction on these 1-D basis functions to obtain a basis of size `(k 2^n)^D`. For a default full grid, this vector space is equivalent to subdividing each axis of this domain into 2^n sub-intervals, for a total of `2^(n D)` hypercubes, and interpolating the function on each of these subdomains by multivariate polynomials of degree strictly less than `k`. The basis functions `V(level, cell, mode)` now have each of their arguments as a `D`-vector.

Because this scheme becomes rapidly and prohibitively expensive in high dimensions, we apply the well-known sparse grid cutoff to exclude all basis functions with levels having 1-norm greater than `n`. This reduces us from `O(k^D 2^(n D))` to `O(k^D 2^n n^(D-1))`.

### Interpolation

The following method interpolates a given function `f` on the `D`-dimensional hypercube `[0,1]^D`.

```julia
coeffs_DG(D::Int, k::Int, n::Int, f::Function; scheme="sparse")
```

For example, to interpolate a 2-D function at order k = 3 and n = 5, we would do

```julia
D = 2; k = 3; n = 5;
f = x->sin(2*pi*x[1])*sin(2*pi*x[2])
full_coeffs   = coeffs_DG(D, k, n, f; scheme="full"  )
sparse_coeffs = coeffs_DG(D, k, n, f; scheme="sparse")
```

This gives a dictionary indexed by multi-levels of type CartesianIndex{2}. For a given multi-level, any dictionary entry is an array of arrays indexed by the multi-cell and multi-mode that holds the coefficient data.

In either case, the interpolation can be evaluated at a given point `xs` (given by an array of length D) by

```julia
full_val   = reconstruct_DG(full_coeffs  , xs)
sparse_val = reconstruct_DG(sparse_coeffs, xs)
```

No `scheme`, `D`, `k`, or `n`, need be specified, as all of these can be deduced from the dictionary itself.

The default scheme is always `"sparse"`.

### Solving Hyperbolic PDEs

All solvers can be found in the `PDEs.jl` script. The simplest one solves the wave equation in 1-D using the DG position basis
	
```julia
wave_evolve_1D(k::Int, n::Int,
			   f0::Function, v0::Function,
			   time0::Real, time1::Real;
			   base = "hier", order = "45").
```
	
Here, `k` and `n` are as before, `f0` is the initial condition for position, `v0` is the initial condition for velocity, `time0` and `time1` are the initial and final times of evolution. For 1D, we can also use the `"pos"` position basis (rather than the multi-resolution), and order `"78"` ode solvers. In either case, `ode45/ode78` from `ODE.jl` are respectively used to solve the differential equation.

For example, for a standing wave solution from `t_0` to `t_1` seconds with k polynomials on each division, up to hierarchical level n, using `ode78`:

```julia
f0 = x->sin(2*pi*x[1])
v0 = x->0
hier_soln = wave_evolve_1D(k, n, f0, v0, t_0, t_1; order="78")
```

For higher dimensions, we use

```julia
wave_evolve(D::Int, k::Int, n::Int,
			f0::Function, v0::Function,
			time0::Real, time1::Real;
			order="45", scheme="sparse")
```

For example, to solve a standing wave in 2-D from `t_0` to `t_1` with polynomials of degree less than `k` on each division, sparse depth `n`, and using `ode78`:

```julia
f0 = x->sin(2*pi*x[1])*sin(2*pi*x[2])
v0 = 0
sparse_soln = wave_evolve(2, k, n, f0, v0, t_0, t_1; order="78")
```

At `n = 1`, this is the same as `wave_evolve_1D` in the `"hier"` basis.

Although the option `scheme = "full"` is implemented for higher dimensions, it is almost never computationally feasible, even in dimensions as low as `D = 3`.

### Coefficient Vectors

It is best to work directly with vectors of coefficients rather than dictionary or tree-like structures when defining linear operators (like differentiation) acting on the interpolations. This way, derivatives are represented by sparse matrices that can act on the coefficient vectors with very fast performance. This is ideal for time evolution.

For this reason we have the `vcoeffs_DG` method that works exactly the same as their non-vector relatives, but return a coefficient vector instead of a dictionary.

We can convert between coefficient vectors and corresponding dictionaries as follows (in 2D with `k = 3, n = 5`):

```julia
D = 2; k = 3; n = 5;
f = x->sin(4*x[1]+x[2])

full_vcoeffs = vcoeffs_DG(D, k, n, f; scheme="full")
full_dict = V2D(D, k, n, full_vcoeffs; scheme="full")

Note this also means full_vcoeffs = D2V(D, k, n, full_dict; scheme="full")

sparse_vcoeffs = vcoeffs_DG(D, k, n, f; scheme="sparse")
sparse_dict = V2D(D, k, n, sparse_vcoeffs; scheme="sparse")

Note this also means sparse_vcoeffs = D2V(D, k, n, sparse_dict; scheme="sparse")
```

In addition, we can generate a lookup by using the `D2Vref` and `V2Dref`that takes us from a dictionary-like level-cell-mode scheme to a specific index in a coefficient vector and vice versa. This is useful when trying to understand specific values in a vector of coefficients.

If we have the vector `sparse_vcoeffs` as above and wanted to understand what `(level, cell, mode)` corresponds to the value of `sparse_vcoeffs[10]`, the code would be

```julia
VD = V2Dref(D, k, n; scheme="sparse")
VD[10]
```

Conversely, If we wanted to see what index we should look at to get the coefficient for `(level, cell, mode) = ((4,4), (1,2), (2, 3))`, the code would be

```julia
level	= CartesianIndex{2}((4, 4))
cell	= CartesianIndex{2}((1, 2))
mode	= CartesianIndex{2}((2, 3))
DV  	= D2Vref(D, k, n; scheme="sparse")
index	= DV[(level, cell, mode)]
```

### Differentiation

The higher dimensional full and sparse derivative operators are implemented in `Multidim_Derivative.jl`. Here, `D`, `k`, and `n` are as before, and `i` specifies the axis along which the derivative is taken.

```julia
D_matrix(D::Int, i::Int, k::Int, n::Int; scheme="sparse")
```
	
We also have the gradient vector that is `D_matrix` over all `i`, and the Laplacian:

```julia
grad_matrix(D::Int, k::Int, n::Int; scheme="sparse")
laplacian_matrix(D::Int, k::Int, n::Int; scheme="sparse")
```

All these matrices act directly on the coefficient vectors obtained from `vcoeffs_DG`. For now, we assume periodic boundary conditions. It is not too difficult to generalize away from a periodic boundary.

In the discontinuous Galerkin method, it is customary to formulate the derivative in the weak sense, and that is how this method is implemented. As is standard, the derivative matrix elements are calculated through an integration by parts, leading to the derivative matrix having two summands: one being a boundary term and the other a derivative term on the interior.

<!-- ## Other Bases

### Canonical "Hat" Basis

To interpolate a multivariate function in the multiresolution (aka hierarchical) scheme of hat-functions (c.f. Gerstner & Griebel below) use:

```julia
coeffs_hat(D::Int, n::Int, f::Function; scheme = "sparse")
```

The coefficients are given as a dictionary of `CartesianIndex{D}` types. The keys of the dictionary correspond to multi-levels and corresponding `CartesianIndex{D}` entries correspond to the set of multi-cells at a specific multi-level.

For example in 2-D:

```julia
coeffs = coeffs_hat(2, n, x->sin(4*x[1]+x[2]); scheme="sparse")
```

For either the full or sparse basis coefficients, we can reconstruct the interpolation at a point xs using the `reconstruct_hat(coefficients::Dict{CartesianIndex{D},Array{Float64,D}}, xs::NTuple{D,T})` function. For example we could take the prior result and evaluate at the point `(.4,.6)` by

```julia
reconstruct_hat(coeffs, (.4,.6))
``` -->

## Construction of the DG Basis

The DG position basis in consists of `k` Legendre polynomials (from degree 0 to `k-1`) supported on subintervals of the form `(i 2^-level, (i+1) 2^-level)`. By construction, this gives an orthonormal basis. For this standard 1-D DG basis, the method `get_vcoeffs(k::Int, level::Int, f::Function)` gives an appropriate array of coefficients that can be reconstructed at a point `x` using `reconstruct_vcoeffs(k,level,coefficients,x)`.

The corresponding multiresolution basis formed by a series of orthogonalizations implemented in `DG_Functions.jl`. This basis spans the same space as the position basis, but makes use of a discontinuous basis (c.f. `1D_DG_Functions.jl`) to achieve the hierarchical structure.

(More to come on orthogonalization schemes)

## Future updates

The family of PDEs that can be solved is currently being extended to include a wide variety of Einstein-like equations. The set of boundary conditions is currently being expanded beyond just periodic and Dirichlet conditions, to account for vanishing conditions at conformal infinity, as well as to handle singularities.

Future updates include adaptivity as well as parallelizability to these methods.

## References

[A Sparse Grid Discontinuous Galerkin Method for high-dimensional transport equations](http://arxiv.org/abs/1602.02124)

[Sparse Grid DG Methods for Elliptic Equations](http://arxiv.org/abs/1508.07781)

Gerstner & Griebel:
[A survey of sparse grid discretization and applciations to quadrature, data compression, and dynamical systems](http://wissrech.ins.uni-bonn.de/research/pub/gerstner/sparsegrids.pdf)

[Introduction to the Interior Penality Discontinuous Galerkin (IPDG)](http://www.cs.elte.hu/~izsakf/otka/presentation_dg_23dim_en.pdf)
