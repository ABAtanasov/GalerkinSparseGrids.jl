# The exact coefficients for f(x_1, ..., x_D) = 1 in the modal basis
function get_one_modal(D::Int, k::Int, n::Int)
	len = get_size(D, k, n)
	return [i==1 ? 1.0 : 0.0 for i in 1:len]
end

# The exact coefficients of the function f(x_1, ..., x_D) = x_i
function get_xi_modal(D::Int, i::Int, k::Int, n::Int)
	len = get_size(1, k, n)
	coeffs1_1D::Array{Float64, 1} = [i==1 ? 1 : 0 for i in 1:len] # Coeffs for 1
	coeffs2_1D::Array{Float64, 1} = [i==2 ? 1/(sqrt(3)) : 0 for i in 1:len]
	coeff_array = [d==i ? coeffs2_1D : coeffs1_1D for d in 1:D]
	return tensor_construct(D, k, n, coeff_array)
end

function get_xi2_modal(D::Int, i::Int, k::Int, n::Int)
	len = get_size(1, k, n)
	# Coefficients for 1:
	coeffs1_1D::Array{Float64, 1} = [i==1 ? 1 : 0 for i in 1:len]
	# Coefficients for (v-1/2)^2:
	coeffs2_1D::Array{Float64, 1} = [i==3 ? 2/(3*sqrt(5)) : 0 for i in 1:len]	 
	coeffs2_1D += 1/3 * coeffs1_1D
	coeff_array = [d==i ? coeffs2_1D : coeffs1_1D for d in 1:D]
	return tensor_construct(D, k, n, coeff_array)
end

# The exact coefficients for f(x) = |x|^2 in the modal basis
function get_r2_modal(D::Int, k::Int, n::Int)
	coeff_array = [get_xi2_modal(D, i, k, n) for i in 1:D]

	return sum(coeff_array)
end

# The exact coefficients for f(x_1, ..., x_D) = x_i = x_i in the point basis
# ^* insofar as the transformation matrices are exact
function get_xi_point(D::Int, i::Int, k::Int, n::Int, 
						m2n::SparseMatrixCSC{Float64, Int},
						n2p::SparseMatrixCSC{Float64, Int})
	return n2p * (m2n * get_xi_modal(D, i, k, n))
end
    
# The exact coefficients for f(x_1, ..., x_D) = |x|^2 in the point basis
# ^* insofar as the transformation matrices are exact
function get_r2_point(D::Int, k::Int, n::Int, 
						m2n::SparseMatrixCSC{Float64, Int},
						n2p::SparseMatrixCSC{Float64, Int})
	v = n2p * (m2n * get_r2_modal(D, k, n))
	v[v .< 0] .= 0
	return v
end