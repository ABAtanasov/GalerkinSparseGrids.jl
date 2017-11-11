#--------------------------------------------------
#
# An explicit tensor constructor of the coefficients
# for \prod_i f_i (x_i) given the coefficients of
# each f(x_i). 
#
#--------------------------------------------------

# This is a useful way to reduce a difficult O(N log^(d-1) N)
# space of integrations down to O(N d) 
# for specifying coefficients of initial conditions

# Efficiency criticality: MEDIUM

# Accuracy criticality: LOW
# Most accuracy it dependent on functions 
# called from other scripts

function tensor_construct(D::Int, k::Int, n::Int, coeffArray; scheme="sparse")
	cutoff = get_cutoff(scheme, D, n)
	coeffs = Dict{CartesianIndex{D}, Array{Array{Float64,D},D}}()
	nD_coeffs = zeros(Float64, D)
	modes = ntuple(i-> k, D)
	ls = ntuple(i-> (n+1), D)
	
	for level in CartesianRange(ls)
		cutoff(level) && continue

		ks = ntuple(i -> 1<<pos(level[i]-2), D)  
		level_coeffs = Array{Array{Float64,D}}(ks)
		lvl = ntuple(i -> level[i]-1,D)
		for cell in CartesianRange(ks)
			cell_coeffs=Array{Float64}(modes)
			for mode in CartesianRange(modes)
				nD_coeffs = [(coeffArray[i])[CartesianIndex{1}((level[i],))][cell[i]][mode[i]] for i in 1:D]
				cell_coeffs[mode] = prod(nD_coeffs)
			end
			level_coeffs[cell]=cell_coeffs
		end
		coeffs[level] = level_coeffs
	end
	return coeffs
end
