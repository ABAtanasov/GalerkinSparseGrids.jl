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
# basis functions in the full or sparse grid schemes
# using the degree k-1 Galerkin polynomials
# with multi-resolution up to n
function get_size(D::Int, k::Int, n::Int; scheme="sparse")
	cutoff = get_cutoff(scheme, D, n)
	ls = ntuple(i->(n+1), D)
	size = 0
	for level in CartesianRange(ls)
		cutoff(level) && continue

		ks = ntuple(q -> 1<<pos(level[q]-2), D)
		size += prod(ks)*k^D
	end
	return size
end


function D2V{d,T<:Real}(D::Int, k::Int, n::Int,
						coeffs::Dict{CartesianIndex{d}, Array{Array{T,d},d}};
						scheme="sparse")
	(d != D) && throw(TypeError(:coeffs))
	cutoff		= get_cutoff(scheme, D, n)
	size		= get_size(D, k, n; scheme=scheme)
	vect		= Array(Float64, size)
	f_numbers	= ntuple(i-> k, D)
	ls			= ntuple(i->(n+1), D)
	j = 1
	for level in CartesianRange(ls) #This really goes from 0 to l_i for each i
		cutoff(level) && continue

		ks = ntuple(q -> 1<<pos(level[q]-2), D)  #This sets up a specific k+1 vector
		for place in CartesianRange(ks)
			for f_number in CartesianRange(f_numbers)
				vect[j] = coeffs[level][place][f_number]
				j+=1
			end
		end
	end
	return vect
end


function V2D{T<:Real}(D::Int, k::Int, n::Int, vect::Array{T}; scheme="sparse")
	cutoff		= get_cutoff(scheme, D, n)
	coeffs		= Dict{CartesianIndex{D}, Array{Array{Float64,D},D}}()
	f_numbers	= ntuple(q-> k, D)
	ls			= ntuple(i->(n+1), D)
	j = 1
	for level in CartesianRange(ls) #This really goes from 0 to l_i for each i
		cutoff(level) && continue

		ks = ntuple(q -> 1<<pos(level[q]-2), D)  #This sets up a specific k+1 vector
		level_coeffs = Array(Array{Float64},ks)	 #all the coefficients at this level
		for place in CartesianRange(ks)
			place_coeffs = Array(Float64,f_numbers)
			for f_number in CartesianRange(f_numbers)
				place_coeffs[f_number] = vect[j]
				j+=1
			end
			level_coeffs[place] = place_coeffs
		end
		coeffs[level] = level_coeffs
	end
	return coeffs
end

function D2Vref(D::Int, k::Int, n::Int; scheme="sparse")
	cutoff		= get_cutoff(scheme, D, n)
	size		= get_size(D, k, n; scheme=scheme)
	dict		= Dict{NTuple{3,CartesianIndex{D}}, Int}()
	f_numbers	= ntuple(q-> k, D)
	ls			= ntuple(i->(n+1), D)
	j = 1
	for level in CartesianRange(ls)
		cutoff(level) && continue

		ks = ntuple(q -> 1<<pos(level[q]-2), D)  #This sets up a specific k+1 vector
		lvl = ntuple(i -> level[i]-1,D)
		for place in CartesianRange(ks)
			for f_number in CartesianRange(f_numbers)
				dict[(level, place, f_number)] = j
				j+=1
			end
		end
	end
	return dict
end

function V2Dref(D::Int, k::Int, n::Int; scheme = "sparse")
	cutoff		= get_cutoff(scheme, D, n)
	size		= get_size(D, k, n; scheme=scheme)
	vect		= Array(NTuple{3,CartesianIndex{D}}, size)
	f_numbers	= ntuple(q-> k, D)
	ls			= ntuple(i->(n+1), D)
	j = 1	
	for level in CartesianRange(ls)
		cutoff(level) && continue

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

function vcoeffs_DG(D::Int, k::Int, n::Int, f::Function;
								rel_tol = REL_TOL, abs_tol = ABS_TOL,
								max_evals=MAX_EVALS,
								scheme="sparse")
	cutoff		= get_cutoff(scheme, D, n)
	len			= get_size(D, k, n; scheme=scheme)
	coeffs		= Array(Float64, len)
	f_numbers	= ntuple(i-> k, D)
	ls			= ntuple(i-> (n+1), D)
	j = 1
	for level in CartesianRange(ls)     # This really goes from 0 to l_i for each i,
		cutoff(level) && continue

		ks = ntuple(i -> 1<<pos(level[i]-2), D)  #This sets up a specific k+1 vector
		lvl = ntuple(i -> level[i]-1,D)
		for place in CartesianRange(ks)
			for f_number in CartesianRange(f_numbers)
				coeffs[j] = get_coefficient_DG(k, lvl, place, f_number, f;
												rel_tol=rel_tol, abs_tol=abs_tol, 
												max_evals=max_evals)
				j+=1
			end
		end
	end
	return coeffs
end

