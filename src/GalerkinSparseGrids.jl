module GalerkinSparseGrids

# package code goes here

include("Hat_Methods.jl")  # Using non-galerkin elementary 'hat' basis functions
include("DG_Functions.jl") # Gram-Schmidt procedure for DG basis functions in 1-D
include("Specific_DG_Functions.jl") # Explicitly building the 1-D Basis
include("DG_Methods.jl") # Multidimensional hierarchical & sparse coefficients
include("DG_Derivative.jl") # 1-D symbolic piecewise derivative 
include("DG_vMethods.jl") # Going between a dictionary & a vector of coeffs
include("DG_Derivative_Matrix.jl") # Precomputing derivative matrix for coeff vect
include("Position_Basis_DG.jl") # Constructing ideal derivative matrix in position space
include("PDEs.jl") # Solving the 1-D wave equation with periodic inital conditions 


#make naming scheme more systematic 
export standard_coefficients
export standard_reconstruct
export hier_coefficients
export sparse_coefficients
export reconstruct 

export hier_coefficients_DG
export sparse_coefficients_DG
export reconstruct_DG

export Full_V2D
export Sparse_V2D
export Full_D2V
export Sparse_D2V
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
export hier_wave_equation
export norm_squared
export pos_energy_func
export hier_energy_func

end # module
