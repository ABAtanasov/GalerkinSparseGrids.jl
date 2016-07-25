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

export total_value



function total_value(n::Int)
	value = 0
	for i in 1:n
		value += 1
	end
	value
end

end # module
