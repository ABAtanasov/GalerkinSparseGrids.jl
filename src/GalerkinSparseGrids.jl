module GalerkinSparseGrids

using HCubature
using StaticArrays
using ODE
using SparseArrays

# The following script files are used:

include("schemes.jl")                    # CartesianIndex Manipulation Methods
include("additional_tools.jl")            # Helper functions not specific to any part of the package
include("dg_basis.jl")                    # Gram-Schmidt procedure for DG basis functions in 1-D
include("1d_dg_functions.jl")             # Explicitly building the 1-D Galerkin Basis
include("dg_methods.jl")                # Multidimensional hierarchical & sparse coefficients
include("dg_vmethods.jl")                # Going between a dictionary & a vector of coeffs
include("derivative_matrix_elements.jl")# 1-D symbolic piecewise derivative
include("derivative_precompute.jl")        # Precomputing derivative matrix for coeff vect
include("1d_derivative.jl")                # Constructing ideal 1D derivative matrix using boundary terms
include("multidim_derivative.jl")        # Multidimensional DG Derivatives in full & sparse bases
include("pdes.jl")                        # Solving the n-D wave equation with periodic boundary
include("1d_nodal_basis.jl")            # Construct 1D nodal basis for multiplication
include("multidim_nodal_basis.jl")        # Construct Multidimensional nodal basis for multiplication

include("basic_function_exact_coeffs.jl") # Exact coefficients for certain easy functions
include("error_measure.jl")                # Monte Carlo Methods to measure error
include("tensor_construct.jl")            # Quickly calculates coeffs of simple tensors of functions

export

coeffs_DG,
reconstruct_DG,

get_size,
V2D,
D2V,
V2Dref,
D2Vref,
vcoeffs_DG,

transform,
make_modal2point_matrices,
make_point2modal_matrices,

D_matrix,
grad_matrix,
laplacian_matrix,

nodal2points_1D,
points2nodal_1D,
nodal2heir_1D,
hier2points_1D,
points2hier_1D,

get_one_modal,
get_xi_modal,
get_r2_modal,
get_xi_point,
get_r2_point,

mcerr,
mcerr2,

tensor_construct,

wave_evolve_1D,
wave_evolve,
vlasov_evolve


end # module
