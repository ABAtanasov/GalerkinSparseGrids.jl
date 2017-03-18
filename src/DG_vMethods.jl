#-------------------------------------------------------------------
# It is advantageous to have the coeffs as
# one single vector in order to make use of BLAS
# and related libraries when defining operators on 
# our space of functions
#
# In this script, we "vectorize" our functions,
# not in the regular meaning of the word, but in the
# sense that dictionaries will all be replaced by vectors
#
# The main difficulty with this is that, unlike with dictionaries
# there is no easy way to go from (level, place, fnumber) to a
# corresponding index in a 1-D vector. 
#
# For this, we have the reference functions:
# full/sparse_referenceV2D/D2V 
# D2V generates a dict that, upon input of a level, place, f_number
# gives the corresponding index in a vector
# V2D generates a vector, with row i having the three numbers
# level, place, f_number corresponding to index i in the vector
# 
# These methods work in all dimensions, and there are ones for both
# full and sparse grids. 
#
# The entire script culminates in a final result: a matrix 
# representation of the derivative operator 
# (both in full and sparse bases)
#-------------------------------------------------------------------


# Calculates the exact dimension of interpolating
# basis functions in the full grid scheme
# using the degree k-1 Galerkin polynomials
# with multi-resolution up to ls
function full_size(k, ls)
	D = length(ls)
	size=0
	f_numbers= ntuple(q -> k, D)
	
	for level in CartesianRange(ls)
		places = ntuple(q -> 1<<pos(level[q]-2), D)
		size += prod(places)*k^D
	end
	return size
end

# Calculates this same dimension when using
# the sparse grid scheme
function sparse_size(k, n, D)
    size=0
    ls = ntuple(i-> (n+1),D)
	
    for level in CartesianRange(ls) 
        diag_level=0;
        for q in 1:D
            diag_level+=level[q]
        end
        if diag_level > n + D 		# If we're past the relevant levels
            continue				# Don't calculate anything
        end  
        ks = ntuple(q -> 1<<pos(level[q]-2), D)
        size += prod(ks)*k^D
    end
    return size
end



function full_D2V{D,T<:Real}(k::Int, coefficients::Dict{CartesianIndex{D}, Array{Array{T,D},D}}, ls::NTuple{D,Int})
	j=1
	size = full_size(k, ls)
	f_numbers= ntuple(q -> k, D)
	vect = Array(Float64,size)

    for level in CartesianRange(ls)     # This really goes from 0 to l_i for each i,
        ks = ntuple(q -> 1<<pos(level[q]-2), D)  #This sets up a specific k+1 vector
        for place in CartesianRange(ks)
			for f_number in CartesianRange(f_numbers)
                vect[j]=coefficients[level][place][f_number]
				j+=1
            end
        end
    end
	return vect
end

function sparse_D2V{D,T<:Real}(k::Int, coefficients::Dict{CartesianIndex{D}, Array{Array{T,D},D}}, n::Int)
	j=1
	size = sparse_size(k,n,D)
	f_numbers= ntuple(i-> k, D)
    ls = ntuple(i->(n+1),D)
	vect = Array(Float64,size)
	
    for level in CartesianRange(ls) #This really goes from 0 to l_i for each i
        diag_level=0;
        for q in 1:D
            diag_level+=level[q]
        end
        if diag_level > n + D #If we're past the levels we care about, don't compute coeffs
            continue
        end  #Otherwise we'll go ahead and DO IT. The same code follows as before.
	    ks = ntuple(q -> 1<<pos(level[q]-2), D)  #This sets up a specific k+1 vector
	    for place in CartesianRange(ks)
            for f_number in CartesianRange(f_numbers)
                vect[j] = coefficients[level][place][f_number]
				j+=1
            end
        end
    end
    return vect
end



function full_V2D{D,T<:Real}(k::Int, vect::Array{T}, ls::NTuple{D,Int})
    coeffs = Dict{CartesianIndex{D}, Array{Array{Float64,D},D}}()
	f_numbers= ntuple(q-> k, D)
	j=1
	
	for level in CartesianRange(ls)     # This really goes from 0 to l_i for each i,
        ks = ntuple(q -> 1<<pos(level[q]-2), D)  #This sets up a specific k+1 vector
        level_coeffs = Array(Array{Float64},ks)	 #all the coefficients at this level
        for place in CartesianRange(ks)
            place_coeffs=Array(Float64,f_numbers)
			for f_number in CartesianRange(f_numbers)
                place_coeffs[f_number]=vect[j]
				j+=1
            end
			level_coeffs[place] = place_coeffs
        end
		coeffs[level] = level_coeffs
    end
	return coeffs
end

function sparse_V2D{T<:Real}(k::Int, vect::Array{T}, n::Int, D::Int)
    coeffs = Dict{CartesianIndex{D}, Array{Array{Float64,D},D}}()
	f_numbers= ntuple(q-> k, D)
    ls = ntuple(i->(n+1),D)
	j=1
	
	for level in CartesianRange(ls) #This really goes from 0 to l_i for each i
        diag_level=0;
        for q in 1:D
            diag_level+=level[q]
        end
        if diag_level > n + D #If we're past the levels we care about, don't compute coeffs
            continue
        end  
		#Otherwise we'll go ahead and DO IT. The same code follows as before.
        ks = ntuple(q -> 1<<pos(level[q]-2), D)  #This sets up a specific k+1 vector
        level_coeffs = Array(Array{Float64},ks)	 #all the coefficients at this level
        for place in CartesianRange(ks)
            place_coeffs=Array(Float64,f_numbers)
			for f_number in CartesianRange(f_numbers)
                place_coeffs[f_number]=vect[j]
				j+=1
            end
			level_coeffs[place]=place_coeffs
        end
		coeffs[level] = level_coeffs
    end
	return coeffs
end



function full_referenceD2V{D}(k::Int, ls::NTuple{D,Int})
	j=1
	size=0
	f_numbers= ntuple(q-> k, D)
	dict = Dict{NTuple{3,CartesianIndex{D}}, Int}()
	
    for level in CartesianRange(ls)
        ks = ntuple(q -> 1<<pos(level[q]-2), D)  #This sets up a specific k+1 vector
		lvl = ntuple(i -> level[i]-1,D)
        for place in CartesianRange(ks)
			for f_number in CartesianRange(f_numbers)
                dict[(level,place,f_number)] = j
				j+=1
            end
        end
    end
	return dict
end

function full_referenceV2D{D}(k::Int, ls::NTuple{D,Int})
	j=1
	size=full_size(k, ls)
	f_numbers= ntuple(q-> k, D)
	vect = Array(NTuple{3,CartesianIndex{D}}, size)
	
    for level in CartesianRange(ls)
        ks = ntuple(q -> 1<<pos(level[q]-2), D)  #This sets up a specific k+1 vector
		lvl = ntuple(i -> level[i]-1,D)
        for place in CartesianRange(ks)
			for f_number in CartesianRange(f_numbers)
                vect[j] = (level, place, f_number)
				j+=1
            end
        end
    end
	return vect
end

function sparse_referenceD2V(k::Int,n::Int,D::Int)
	j=1
	f_numbers= ntuple(q-> k, D)
	size=sparse_size(k,n,D)
	ls = ntuple(i->(n+1),D)
	dict = Dict{NTuple{3,CartesianIndex{D}}, Int}()
	
    for level in CartesianRange(ls)
        diag_level=0;
        for q in 1:D
            diag_level+=level[q]
        end
        if diag_level > n + D #If we're past the levels we care about, don't compute coeffs
            continue
		end
        ks = ntuple(q -> 1<<pos(level[q]-2), D)  #This sets up a specific k+1 vector
		lvl = ntuple(i -> level[i]-1,D)
        for place in CartesianRange(ks)
			for f_number in CartesianRange(f_numbers)
                dict[(level,place,f_number)] = j
				j+=1
            end
        end
    end
	return dict
end

function sparse_referenceV2D(k::Int,n::Int,D::Int)
	f_numbers= ntuple(q-> k, D)
	size=sparse_size(k,n,D)
	vect = Array(NTuple{3,CartesianIndex{D}}, size)
    ls = ntuple(i->(n+1),D)
	j=1
	
    for level in CartesianRange(ls)
        diag_level=0;
        for q in 1:D
            diag_level+=level[q]
        end
        if diag_level > n + D #If we're past the levels we care about, don't compute coeffs
            continue
		end
        ks = ntuple(q -> 1<<pos(level[q]-2), D)  #This sets up a specific k+1 vector
		lvl = ntuple(i -> level[i]-1,D)
        for place in CartesianRange(ks)
			for f_number in CartesianRange(f_numbers)
                vect[j] = (level, place, f_number)
				j+=1
            end
        end
    end
	return vect
end


#------------------------------------------------------
# Let's now make the coefficient operators work on
# and return vectors
#------------------------------------------------------

function vhier_coefficients_DG{D}(k::Int, f::Function, ls::NTuple{D,Int};
							rel_tol = REL_TOL, abs_tol = ABS_TOL, max_evals=MAX_EVALS)
 	l = k^D * 2^(sum(ls)-D)
    coeffs = Array(Float64, l)
	f_numbers= ntuple(i-> k, D)
	j=1
    for level in CartesianRange(ls)     # This really goes from 0 to l_i for each i,
        ks = ntuple(i -> 1<<pos(level[i]-2), D)  #This sets up a specific k+1 vector
		lvl = ntuple(i -> level[i]-1,D)
        for place in CartesianRange(ks)
			for f_number in CartesianRange(f_numbers)
                coeffs[j]=get_coefficient_DG(k,f,lvl,place,f_number;
									rel_tol = rel_tol, abs_tol=abs_tol, max_evals=max_evals)
				j+=1          
            end
        end
    end
    return coeffs
end



function vsparse_coefficients_DG(k::Int, f::Function, n::Int, D::Int;
								rel_tol = REL_TOL, abs_tol = ABS_TOL, max_evals=MAX_EVALS)
 	len = sparse_size(k,n,D)
    coeffs = Array(Float64, len)
	f_numbers= ntuple(i-> k, D)
	ls = ntuple(i-> (n+1),D)
	j=1
    for level in CartesianRange(ls)     # This really goes from 0 to l_i for each i,
        diag_level=0;
        for q in 1:D
            diag_level+=level[q]
        end
        if diag_level > n + D #If we're past the levels we care about, don't compute coeffs
            continue
        end  
        ks = ntuple(i -> 1<<pos(level[i]-2), D)  #This sets up a specific k+1 vector
		lvl = ntuple(i -> level[i]-1,D)
        for place in CartesianRange(ks)
			for f_number in CartesianRange(f_numbers)
                coeffs[j]=get_coefficient_DG(k,f,lvl,place,f_number;
									rel_tol = rel_tol, abs_tol=abs_tol, max_evals=max_evals)
				j+=1          
            end
        end
    end
    return coeffs
end

# function full_eval{D}(k::Int, ls::NTuple{D,Int})
# 	j = 1
# 	size = full_size(k, ls)
# 	f_numbers = ntuple(q-> k, D)
# 	vect = zeros(T, size)
# 	lookup = full_referenceV2D(k, ls)
# 	dict = full_V2D(k, vect, ls)		# dict of 0.0s
#
# 	for i in 1:size
# 		level, place, number = lookup[i]
# 		dict[level][place][number] = 1.0
#
# 		dict[level][place][number] = 0.0
# 	end
# end
#
# function sparse_eval{D}(k::Int, vect::Array{T}, ls::NTuple{D,Int})
#
# end

# You may want to implement fullreconstruct and sparsereconstruct

