module GalerkinSparseGrids

# The prerequisite packages as of May 2018 are Cubature.jl and ODE.jl

using HCubature
using StaticArrays
using ODE
using SparseArrays

# The following script files are used

include("Schemes.jl")					# CartesianIndex Manipulation Methods
include("Additional_Tools.jl")			# Helper functions not specific to any part of the package
# include("Hat_Methods.jl")				# Using non-galerkin elementary 'hat' basis functions
include("DG_Basis.jl")					# Gram-Schmidt procedure for DG basis functions in 1-D
include("1D_DG_Functions.jl") 			# Explicitly building the 1-D Galerkin Basis
include("DG_Methods.jl")				# Multidimensional hierarchical & sparse coefficients
include("DG_vMethods.jl")				# Going between a dictionary & a vector of coeffs
include("Derivative_Matrix_Elements.jl")# 1-D symbolic piecewise derivative
include("Derivative_Precompute.jl")		# Precomputing derivative matrix for coeff vect
include("1D_Derivative.jl")				# Constructing ideal 1D derivative matrix using boundary terms
include("Multidim_Derivative.jl")		# Multidimensional DG Derivatives in full & sparse bases
include("PDEs.jl")						# Solving the n-D wave equation with periodic boundary
include("1D_Nodal_Basis.jl")			# Construct 1D nodal basis for multiplication
include("Multidim_Nodal_Basis.jl")		# Construct Multidimensional nodal basis for multiplication

include("Error_Measure.jl")				# Monte Carlo Methods to measure error
include("Tensor_Construct.jl")			# Quickly calculates coeffs of simple tensors of functions
include("Traveling_Wave_Example.jl")	# Construct and evolve traveling waves

export

coeffs_DG,
reconstruct_DG,

get_size,
V2D,
D2V,
V2Dref,
D2Vref,
vcoeffs_DG,

transform_1D,

D_matrix,
grad_matrix,
laplacian_matrix,

wave_evolve_1D,
wave_evolve,

nodal2points_1D,
points2nodal_1D,
nodal2heir_1D,
hier2points_1D,
points2hier_1D,

mcerr,
mcerr2,

tensor_construct,

cos_coeffs,
sin_coeffs,
traveling_wave_solver


end # module
