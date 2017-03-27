module GalerkinSparseGrids

# The prerequisite packages as of January 2017 are Cubature.jl and ODE.jl

using Cubature
using ODE

# The following script files are used

include("Schemes.jl")				   # CartesianIndex Manipulation Methods
include("Hat_Methods.jl") 			   # Using non-galerkin elementary 'hat' basis functions
include("DG_Basis.jl")	  			   # Gram-Schmidt procedure for DG basis functions in 1-D
include("1D_DG_Functions.jl") 		   # Explicitly building the 1-D Galerkin Basis
include("DG_Methods.jl") 			   # Multidimensional hierarchical & sparse coefficients
include("DG_vMethods.jl") 			   # Going between a dictionary & a vector of coeffs
include("DG_Derivative_Matrix_Elements.jl") 		   # 1-D symbolic piecewise derivative 
include("DG_Derivative_Precompute.jl") # Precomputing derivative matrix for coeff vect
include("Derivative_LF_1D.jl") 	   	   # Constructing ideal 1D derivative matrix using boundary terms
include("Multidim_Derivative.jl") 	   # Multidimensional DG Derivatives in full & sparse bases
include("PDEs.jl") 					   # Solving the n-D wave equation with periodic boundary


include("Error_Measure.jl") 		   # Monte Carlo Methods to measure error
include("Tensor_Construct.jl") 		   # Quickly calculates coeffs of `simple tensor' functions
include("Traveling_Wave_Example.jl")


# make naming scheme more systematic 
export standard_coeffs
export standard_reconstruct
export coeffs_hat
export reconstruct_hat 

export coeffs_DG
export reconstruct_DG

export V2D
export D2V
export total_value
export V2Dref
export D2Vref
export vcoeffs_DG

export get_coeffs
export get_vcoeffs
export reconstruct_coeffs
export reconstruct_vcoeffs

export pos_DLF_Matrix
export periodic_pos_DLF_Matrix
export hier_DLF_Matrix
export periodic_hier_DLF_Matrix
export D_matrix
export grad_matrix
export laplacian_matrix

export wave_evolve_1D
export wave_evolve

export mcerr
export mcerr2

export tensor_construct_full
export tensor_construct_sparse

export cos_coeffs
export sin_coeffs
export traveling_wave_solver


end # module
