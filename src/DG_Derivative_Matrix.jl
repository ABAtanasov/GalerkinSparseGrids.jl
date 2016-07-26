#------------------------------------------------------
# Precompute relevant 1-D basis elements
#------------------------------------------------------


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