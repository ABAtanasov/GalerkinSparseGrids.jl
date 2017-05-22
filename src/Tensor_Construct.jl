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
	c = zeros(Float64, D)
	f_numbers= ntuple(i-> k, D)
	ls = ntuple(i-> (n+1), D)
	
	for level in CartesianRange(ls)
		cutoff(level) && continue

		ks = ntuple(i -> 1<<pos(level[i]-2), D)  
		level_coeffs = Array(Array{Float64,D},ks)
		lvl = ntuple(i -> level[i]-1,D)
		for place in CartesianRange(ks)
			place_coeffs=Array(Float64,f_numbers)
			for f_number in CartesianRange(f_numbers)
				c = [(coeffArray[i])[CartesianIndex{1}((level[i],))][place[i]][f_number[i]] for i in 1:D]
				place_coeffs[f_number] = prod(c)
			end
			level_coeffs[place]=place_coeffs
		end
		coeffs[level] = level_coeffs
	end
	return coeffs
end
