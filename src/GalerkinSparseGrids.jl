module GalerkinSparseGrids

# The prerequisite packages as of March 2017 are Cubature.jl and ODE.jl

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

export 
standard_coeffs,
standard_reconstruct,
coeffs_hat,
reconstruct_hat,

coeffs_DG,
reconstruct_DG,

V2D,
D2V,
V2Dref,
D2Vref,
vcoeffs_DG, 

periodic_DLF_Matrix,
D_matrix,
grad_matrix,
laplacian_matrix,

wave_evolve_1D,
wave_evolve,

mcerr,
mcerr2,

tensor_construct,

cos_coeffs,
sin_coeffs,
traveling_wave_solver


end # module
