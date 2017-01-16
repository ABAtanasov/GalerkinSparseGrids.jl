#------------------------------------------------------
# Precompute relevant 1-D basis elements
#------------------------------------------------------

# Efficiency criticality: MED
# Computation done once, but is done upon import
# This could be improved

KMAX = 5
LMAX = 10

precomputed_diffs = Dict{NTuple{4,Int},Array{Float64,2}}()

for level in 0:10
    for place in 1:(1<<pos(level-1)) 
        for f_number in 1:2
            precomputed_diffs[(2,level,place,f_number)] = diff_basis_DG(2,level,place,f_number)
        end
    end
end

for level in 0:10
    for place in 1:(1<<pos(level-1)) 
        for f_number in 1:3
            precomputed_diffs[(3,level,place,f_number)] = diff_basis_DG(3,level,place,f_number)
        end
    end
end

for level in 0:10
    for place in 1:(1<<pos(level-1)) 
        for f_number in 1:4
            precomputed_diffs[(4,level,place,f_number)] = diff_basis_DG(4,level,place,f_number)
        end
    end
end

for level in 0:10
    for place in 1:(1<<pos(level-1)) 
        for f_number in 1:5
            precomputed_diffs[(5,level,place,f_number)] = diff_basis_DG(5,level,place,f_number)
        end
    end
end


# hier2pos = Array(SparseMatrixCSC{Float64,Int64}, (KMAX, LMAX))
#
#
# # NOTE here this precomputation is done with ABS_TOL default.
# # Will want to change in the future
# for k in 1:KMAX
# 	for level in 1:LMAX
# 		hier2pos[k, level] = hier2pos_precompute(k, level)
# 	end
# end
#
#
