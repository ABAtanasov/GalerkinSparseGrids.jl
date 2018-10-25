# -------------------------------------------------
#
# An explicit tensor constructor of the coefficients
# for \prod_i f_i (x_i) given the coefficients of
# each f(x_i).
#
# -------------------------------------------------

# This is a useful way to reduce a difficult O(N log^(d-1) N)
# space of integrations down to O(N d)
# for specifying coefficients of initial conditions

# Efficiency criticality: MEDIUM

# Accuracy criticality: LOW
# Most of the accuracy dependent on functions
# called from other scripts

function tensor_construct(D::Int, k::Int, n::Int, coeff_array::Array{Dict{CartesianIndex{1},
	Array{Array{T, 1}, 1}}, 1}; scheme = "sparse") where T <: Real

	cutoff = get_cutoff(scheme, D, n)
	coeffs = Dict{CartesianIndex{D}, Array{Array{T,D},D}}()
	# nD_coeffs = zeros(T, D)
	ls::NTuple{D, Int}    = ntuple(i-> (n+1), D)
	modes::NTuple{D, Int} = ntuple(i-> k, D)

	for level in CartesianRange(ls)
		cutoff(level) && continue

		cells::NTuple{D, Int} = ntuple(i -> 1<<max(0, level[i]-2), D)
		level_coeffs = Array{Array{T,D}}(cells)
		for cell in CartesianRange(cells)
			cell_coeffs = Array{T}(undef, modes)
			for mode in CartesianRange(modes)
				val = one(T)
				for d in 1:D
					coeff = coeff_array[d]
					tup::NTuple{1, Int} = (level[d],) # This is somehow still slow
					l = CartesianIndex{1}(tup);
					c = cell[d];
					m = mode[d]
					val *= coeff[l][c][m]
				end
				cell_coeffs[mode] = val
			end
			level_coeffs[cell] = cell_coeffs
		end
		coeffs[level] = level_coeffs
	end
	return coeffs
end
