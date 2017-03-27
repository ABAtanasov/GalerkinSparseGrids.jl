#------------------------------------------------------
# Precompute relevant 1-D basis elements
#------------------------------------------------------

# Efficiency criticality: MED
# Computation done once, but is done upon import
# This could be improved

KMAX = 5
LMAX = 10

precomputed_diffs = Dict{NTuple{4,Int},Array{Float64,2}}()

for k in 1:KMAX
	for level in 0:LMAX
	    for place in 1:(1<<pos(level-1)) 
	        for f_number in 1:k
	            precomputed_diffs[(k,level,place,f_number)] = diff_basis_DG(k,level,place,f_number)
	        end
	    end
	end
end
