# -----------------------------------------------------------
#
# Methods for DG interpolation using readable
# dictionary-style structs
#
# -----------------------------------------------------------

using Cubature

const REL_TOL = 1.0e-8
const ABS_TOL = 1.0e-12
const MAX_EVALS = 1000

# -----------------------------------------------------
# Here we have all our methods for
# taking a function and returning an appropriate list
# of DG coefficients.
#
# -----------------------------------------------------

# Efficiency criticality: HIGH

# -----------------------------------------------------
# 1-D Basis
# -----------------------------------------------------

function v{T<:Real}(k::Int, level::Int, cell::Int, mode::Int, x::T)
	if level==0 # At the base level, we are not discontinuous, and we simply
				# use Legendre polynomials up to degree k-1 as a basis
		return LegendreP(mode-1,2*x-1)*sqrt(2.0)
	else
		return h(k, mode, (1<<level)*x - (2*cell-1)) * sqrt(one(T)*(1<<level))
		# Otherwise we use an appropriately shifted, scaled, and normalized
		# DG function
	end
end

function v(k::Int, level::Int, cell::Int, mode::Int)
	return (xs-> v(k,level,cell,mode,xs))
end

# -----------------------------------------------------
# Tensor Product Construction
# -----------------------------------------------------

# Returns the value of the function at x
function V{D,T<:Real}(k::Int, level::NTuple{D,Int},
	cell::CartesianIndex{D}, mode::CartesianIndex{D}, xs::Array{T,1})
	ans = one(T)
	for i = 1:D
		ans *= v(k, level[i], cell[i], mode[i], xs[i])
	end
	return ans
end
# Is there any fast way to precompute the ones I care about??

# Returns a function
function V{D}(k, level::NTuple{D,Int}, 
			cell::CartesianIndex{D}, mode::CartesianIndex{D})
	return (xs-> V(k, level, cell, mode, xs))
end

# -----------------------------------------------------
# Methods for obtaining the coefficients
# -----------------------------------------------------

# Given a 1-D position and level, this tells us which cell 
# that position belongs to, at that level resolution
function cell_index(x::Real,l::Int)
	if l <= 1
		return 1
	end
	if x >= 1
		return 2^(l-1)
	else
		return 1 + floor(Int, 2^(l-1) * x)
	end
end


# This takes an inner product, but since for higher levels our inner product
# is only concerned with a specific region in the grid, we restrict to that
# appropriately, depending on the level
function inner_product{D}(f::Function, g::Function, lvl::NTuple{D,Int}, cell::CartesianIndex{D};
								rel_tol = REL_TOL, abs_tol = ABS_TOL, max_evals=MAX_EVALS)
								
	h = (x-> f(x)*g(x))
	xmin = ntuple(i-> (cell[i]-1)/(1<<max(0,lvl[i]-1)), D)
	xmax = ntuple(i-> (cell[i])/(1<<max(0,lvl[i]-1)), D)
	val = hcubature(h, xmin, xmax; reltol=rel_tol, abstol=abs_tol, maxevals=max_evals)[1]
	return val 
end

# We obtain coefficients simply by doing inner products, it's easy :)
# Only hard part is inner product integrations can be slower than we want :(
function get_coefficient_DG{D}(k::Int, lvl::NTuple{D,Int},
									cell::CartesianIndex{D},
									mode::CartesianIndex{D},
									f::Function;
									rel_tol = REL_TOL,
									abs_tol = ABS_TOL,
									max_evals=MAX_EVALS)

	return inner_product(f, V(k,lvl,cell,mode),lvl,cell; 
							rel_tol = rel_tol, abs_tol = abs_tol, max_evals=max_evals)
end 



# -----------------------------------------------------
# Full or Sparse Galerkin Coefficients in n-D
# -----------------------------------------------------
function coeffs_DG(D::Int, k::Int, n::Int, f::Function; 
										rel_tol = REL_TOL, abs_tol = ABS_TOL,
										max_evals=MAX_EVALS,
										scheme="sparse")
	cutoff	= get_cutoff(scheme, D, n)
	coeffs	= Dict{CartesianIndex{D}, Array{Array{Float64,D},D}}()
	modes	= ntuple(i-> k, D)
	ls		= ntuple(i->(n+1),D)

	for level in CartesianRange(ls) #This really goes from 0 to l_i for each i
		cutoff(level) && continue
		
		cells = ntuple(i -> 1<<max(0, level[i]-2), D)
		level_coeffs = Array{Array{Float64,D}}(cells)
		lvl = ntuple(i -> level[i]-1,D)
		for cell in CartesianRange(cells)
			cell_coeffs = Array{Float64}(modes)
			for mode in CartesianRange(modes)
				cell_coeffs[mode] = get_coefficient_DG(k, lvl, cell, mode, f)
			end
			level_coeffs[cell]=cell_coeffs
		end
		coeffs[level] = level_coeffs
	end
	return coeffs
end


# -----------------------------------------------------------
# Reconstruction (full and sparse) in n-D from a Dict of
# coefficients
# -----------------------------------------------------------

function reconstruct_DG{D,T<:Real}(coeffs::Dict{CartesianIndex{D}, Array{Array{Float64,D},D}},
									xs::Array{T,1})
	
	value	= zero(T)
	k		= size(first(values(coeffs))[1])[1]
	modes	= ntuple(i-> k ,D)

	for key in keys(coeffs)
		level = ntuple(i->key[i]-1,D)
		cell = CartesianIndex{D}(ntuple(i->cell_index(xs[i],level[i]),D))
		coeff = coeffs[key][cell]::Array{T,D}
		@inbounds for mode in CartesianRange(modes)
			value += coeff[mode]*V(k, level, cell, mode, xs)
		end
	end
	return value
end

