#------------------------------------------------------
#
# Sparse Integration Techniques for Calculating Error
#
#------------------------------------------------------

# Efficiency criticality: MEDIUM
# Important for calculating error for testing package
# Not critical for user
# Potential use to replace multidimensional integration 
# from hcubature with more efficient method

function hypervol{D}(ks::NTuple{D,Int})
	ans = one(Float64)
	for i in 1:D
		ans *= ks[i]
	end
	return 1/ans
end

function squadrature(f::Function, n::Int, D::Int)
	coeffs = sparse_coefficients(f, n ,D)
	ans = zero(Float64)
    ls = ntuple(i -> (n+2), D)
	
	for level in CartesianRange(ls)
        diag_level = 0;
        for i in 1:D
            diag_level += level[i]
        end
        if diag_level > n + 2*D 		
            continue					
        else  							
			
            ks = ntuple(i -> 1<<pos(level[i]-3), D)
			block = zeros(Float64)
            for place in CartesianRange(ks)
                block += 0.5*coeffs[level][place] 
            end
			block *= hypervol(ks)
			ans += block
        end
	end
	
	return ans
end