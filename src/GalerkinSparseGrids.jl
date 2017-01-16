module GalerkinSparseGrids

# The prerequisite packages as of January 2017 are Cubature.jl and ODE.jl

using Cubature
using ODE

# The following script files are used

include("Hat_Methods.jl") 			   # Using non-galerkin elementary 'hat' basis functions
include("DG_Basis.jl")	  			   # Gram-Schmidt procedure for DG basis functions in 1-D
include("SQuadrature.jl")  			   # Gram-Schmidt procedure for DG basis functions in 1-D
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


#make naming scheme more systematic 
export standard_coefficients
export standard_reconstruct
export hier_coefficients
export sparse_coefficients
export reconstruct 

export hier_coefficients_DG
export sparse_coefficients_DG
export reconstruct_DG

export full_V2D
export sparse_V2D
export full_D2V
export sparse_D2V
export total_value
export full_referenceD2V
export full_referenceV2D
export sparse_referenceD2V
export sparse_referenceV2D
export vhier_coefficients_DG
export vsparse_coefficients_DG
export sD_matrix

export get_coeffs
export get_vcoeffs
export reconstruct_coeffs
export reconstruct_vcoeffs
export pos_DLF_Matrix
export periodic_pos_DLF_Matrix
export hier_DLF_Matrix
export periodic_hier_DLF_Matrix

export pos_wave_equation4
export pos_wave_equation45
export hier_wave_equation45

export full_D_matrix
export sparse_D_matrix

export sparse_wave_equation45
export sparse_wave_equation78

export squadrature
export mcerr
export mcerr2

export tensor_construct_full
export tensor_construct_sparse

export cos_coeffs
export sin_coeffs
export traveling_wave_equation45
export traveling_wave_equation78


end # module
