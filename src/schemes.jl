# -----------------------------------------------------
# This script gives the methods for maniupulating the
# Cartesian indices corresponding to multilevels that
# are either kept or cut in various schemes of
# "sparsification"
# -----------------------------------------------------

# Efficiency Criticality: LOW
# This is not evaluated very often

# Accuracy Criticality: N/A
# No float manipulation



# Gives the appropriate boolean cutoff corresponding
# to a given scheme, e.g. sparse basis, full basis,
# and in the future possibly the energy basis of 
# Bungartz and Griebel

function get_cutoff(scheme::String, D::Int, n::Int)
	if scheme == "sparse"
		return x -> (sum(x.I) > n+D)
	elseif scheme == "full"
		return x -> false
	else
		throw(ArgumentError)
	end
end