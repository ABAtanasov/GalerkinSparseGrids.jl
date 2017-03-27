#------------------------------------------------------
# This script gives the methods for maniupulating the
# Cartesian indices corresponding to multilevels that
# are either kept or cut in various schemes of
# "sparsification"
#------------------------------------------------------

# Efficiency Criticality: LOW
# This is not evaluated very often

# Accuracy Criticality: N/A
# No float manipulation


# Currently there is an inefficiency in Julia v0.5 with ntuple
# that will hopefully be fixed in the next version. To avoid
# slowdown, we implement a tuple constructor independent of Julia's
# metaprogramming. 
#
# The method below is used for building the multi-dimensional 
# derivative matrix:

@generated function make_cartesian_index{D}(i::Int, arr1::Int, arr2::CartesianIndex{D})
	levs = [:($q == i ? arr1 : arr2[$q]) for q=1:D]
	quote
		$(Expr(:meta, :inline))
		CartesianIndex{D}($(levs...))
	end
end

# Sum method for CartesianIndex class

function csum{D}(level::CartesianIndex{D})
	ans = 0
	for i in 1:D
		ans += level[i]
	end
	return ans
end

# Gives the appropriate boolean cutoff corresponding
# to a given scheme, e.g. sparse basis, full basis,
# and in the future possibly the energy basis of 
# Bungartz and Griebel

function get_cutoff(scheme::String, D::Int, n::Int)
	if scheme == "sparse"
		return x -> (csum(x) > n+D)
	elseif scheme == "full"
		return x -> false
	else
		throw(ArgumentError)
	end
end