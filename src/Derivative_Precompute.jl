# -----------------------------------------------------
# Precompute relevant 1-D basis elements
# -----------------------------------------------------

# Efficiency criticality: MED
# Computation done once, but is done upon import
# This could be improved

KMAX = 5
LMAX = 10

precomputed_diffs = Dict{NTuple{4, Int}, Array{Float64, 2}}()

function precompute_diffs()
	if length(precomputed_diffs) == 0
		print("Precomputing 1D derivative matrix elements... ")
		for k in 1:KMAX
			for level in 0:LMAX
			    for cell in 1:(1<<max(0, level-1)) 
			        for mode in 1:k
			            precomputed_diffs[(k,level,cell,mode)] = diff_basis_DG(k,level,cell,mode)
			        end
			    end
			end
		end
		println("done.")
		return true
	end
	return false
end