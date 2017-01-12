#--------------------------------------------------
#
# An explicit tensor constructor of the coefficients
# for \prod_i f_i (x_i) given the coefficients of
# each f(x_i). 
#
#--------------------------------------------------

# This is a useful way to reduce a difficult O(N log^(d-1) N)
# space of integrations down to O(N d)


# Efficiency criticality: MEDIUM

function tensor_construct_full{D}(k::Int, ls::NTuple{D,Int}, coeffArray)
    coeffs = Dict{CartesianIndex{D}, Array{Array{Float64,D},D}}()
    c = zeros(Float64, D)
	f_numbers= ntuple(i-> k, D)
    for level in CartesianRange(ls)     
        ks = ntuple(i -> 1<<pos(level[i]-2), D)  
        level_coeffs = Array(Array{Float64,D},ks)
		lvl = ntuple(i -> level[i]-1,D)
        for place in CartesianRange(ks)
            place_coeffs=Array(Float64,f_numbers)
			for f_number in CartesianRange(f_numbers)
                c = [coeffArray[i][CartesianIndex{1}((level[i],))][place[i]][f_number[i]] for i in 1:D]
                place_coeffs[f_number]=prod(c)     
            end
            level_coeffs[place]
        end
        coeffs[level] = level_coeffs
    end
    return coeffs
    
end

function tensor_construct_full(k::Int, ls::NTuple{2,Int}, coeffs1)
    return tensor_construct_full(k, ls::NTuple{2,Int}, coeffs1, coeffs1)
end



function tensor_construct_sparse(k::Int, n::Int, D::Int, coeffArray)
    coeffs = Dict{CartesianIndex{D}, Array{Array{Float64,D},D}}()
    c = zeros(Float64, D)
	f_numbers= ntuple(i-> k, D)
    ls = ntuple(i->n+1, D)
    for level in CartesianRange(ls)
        diag_level=0
        for i in 1:D
            diag_level+=level[i]
        end
        if diag_level > n + D  	
            continue			
        end  
		
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

function tensor_construct_sparse(k::Int, n::Int, coeffs1)
    return tensor_construct_sparse(k, n::Int, coeffs1, coeffs1)
end