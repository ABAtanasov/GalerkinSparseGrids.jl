# -----------------------------------------------------------
#
# Evaluating the 1D Galerkin Basis functions
#
# -----------------------------------------------------------

const K_max = 10;

# Efficiency criticality: HIGH
# Central to the efficiency of the package

# Performing fastmath on an array of coefficients
# of x^k and |x|^k, giving the value at a specified point
# supported on [-1, 1]
@fastmath function array2poly(v::Array{T, 1}, x::Real) where T <: Real
	if abs(x)>1
		return zero(T)
	end
	n = length(v)
	k = div(n,2)
	s = zero(T)
	xi = one(T)
	sgn = sign(x)
	@inbounds for i in k:-1:1
		# Using Horner's method
		s *= x
		s += v[i] + flipsign(v[i+k], x)#v[i+k]*sign(x)#
	end
	return s
end

function array2poly(v::Array{T, 1}) where T <: Real
	return (x-> array2poly(v,x))
end

#precompute Legendre
leg_coeffs = legendre(K_max)

# Legendre polynomial supported on [-1, 1]
function LegendreP(k, x)
	k<=K_max || throw(DomainError())
	return array2poly(leg_coeffs[k+1], x)
end

function LegendreP(k)
	k<=K_max || throw(DomainError())
	return array2poly(leg_coeffs[k+1])
end

# precomputing the DG functions
# TODO Make this a 2D array:
#dg_coeffs = Array{Array{Array{Float64,1},1}}(K_max)
dg_coeffs = Array{Vector{Vector{Float64}}}(undef, K_max) # Thanks to Alex Arslan for this

for i in 1:K_max
	dg_coeffs[i] = dg_basis(i)
end

# This is the dg basis function corresponding a given mode
# supported on [-1, 1]
function h(k, mode, x)
	mode<=k || throw(DomainError())
	return array2poly((dg_coeffs[k])[mode],x)
end

function h(k, mode)
	mode<=k || throw(DomainError())
	return array2poly((dg_coeffs[k])[mode])
end

# ----------------------------------------------
#
# Methods relating to the 1D Position Basis:
#
# ----------------------------------------------

# Legendre polynomial on [0, 1]
function leg(mode::Int, x::T) where T <: Real
	return sqrt(2.0)*LegendreP(mode-1, 2*x-1)
end

# level >= 1
# cell in 1:1<<level
# mode in 1:k
# Legendre polynomial on [(cell-1)*h, cell*h] with h = 1/(1<<level)
function basis(level::Int, cell::Int, mode::Int, x::T) where T <: Real
	return leg(mode, (1<<level)*x - (cell-1)) * (2.0)^(level/2)
end

function basis(level::Int, cell::Int, mode::Int)
	return x->basis(level, cell, mode, x)
end

# Perform numerical integration to get position basis coeffs
#
# Never used except if user demands to evolve wave equation
# using only the position basis, never going through heir
function pos_vcoeffs_DG(k::Int, level::Int, f::Function;
						rel_tol = REL_TOL, abs_tol = ABS_TOL, max_evals = MAX_EVALS)
	vcoeffs = Array{Float64}(undef, (1<<level)*(k))
	i = 1
	for cell in 1:(1<<level)
		for mode in 1:k
			fcn			= x->(basis(level, cell, mode,x)*f(x))
			left_bound  = (cell-1)/(1<<level)
			right_bound = (cell)/(1<<level)
			vcoeffs[i]  = hquadrature(fcn, left_bound, right_bound; abstol=abs_tol)[1]
			i += 1
		end
	end
	return vcoeffs
end

# ----------------------------------------------
# Building the hier2pos matrix:
# ----------------------------------------------

# Given an array of the type above, of coefficients
# for both both x^k and |x|^k, then given a choice of
# side (right vs. left) this becomes just a pure
# polynomial array of half the length
function convert_polyarray(v::Array{T, 1}, side = "left") where T <: Real
    n = Int(length(v)//2)
	if side == "left"
		return [v[i]-v[n+i] for i in 1:n]
	elseif side == "right"
		[v[i]+v[n+i] for i in 1:n]
	else throw(ArgumentError())
	end
end

# Gives the polynomial array corresponding to p(x/c)
function scale_polyarray(v::Array{T, 1}, c::Real) where T <: Real
	n = length(v)
	return [v[i] * (1//c^(i-1)) for i in 1:n]
end

# Gives the polynomial array corresponding to p(x-a)
function shift_polyarray(v::Array{T, 1}, a::Real) where T <: Real
	n = length(v)
	v_new = zeros(v)
	for i::Int in 1:n
		for j::Int in 1:i
			v_new[j] += v[i] * (-a)^(i-j) * binomial(i-1, j-1)
		end
	end
	return v_new
end

# Performs symbolic integration of the polynomials
# represented by v and w on the interval [a, b]
function integrate_polyarray(v::Array{T, 1}, w::Array{T, 1}; a::Real = 0, b::Real = 1) where T <: Real
	ans = zero(T)
	for i::Int in 1:length(v)
		for j::Int in 1:length(w)
			ans += v[i] * w[j] * ((b^(i+j-1)-a^(i+j-1))//(i+j-1))
		end
	end
	return ans
end

# With a bit of case work, we can do an explicit
# symbolic integration of the elements in the
# hierarchical basis against the legendre polynomials
# defining the position basis:
function hier2pos(k::Int, max_level::Int, level::Int, cell::Int, mode::Int)

	ans = Real[]

	# Compute the support and discontinuity
	# of the hierarchical basis element
	left_point  = (cell-1)//(1<<max(0, level-1))
	mid_point   = (cell- 1//2)//(1<<max(0, level-1))
	right_point = cell//(1<<max(0, level-1))

	# Level 0 of the hierarchical basis is a continous
	# Legendre polynomial
	if level == 0
		vl = vr = sqrt(2)*shift_polyarray(
								scale_polyarray(leg_coeffs[mode][1:k],
												1//2),
											1//2)

	# The higher DG levels are given by two polynomials
	# defined on neighboring (right & left) intervals of size h/2,
	# which would need to be integrated against individually
	else
		h = 1//(1<<max(0, level-1))
		vl = shift_polyarray(
				scale_polyarray(
					convert_polyarray(dg_coeffs[k][mode], "left"),
								1//2),
							1//2)
		vl = shift_polyarray(scale_polyarray(vl, h), (cell-1)*h)*sqrt(1<<(level))

		vr = shift_polyarray(
				scale_polyarray(
					convert_polyarray(dg_coeffs[k][mode], "right"),
					 			1//2),
							1//2)
		vr = shift_polyarray(scale_polyarray(vr, h), (cell-1)*h)*sqrt(1<<(level))
	end

	h = 1//(1<<max_level)
	for i in 1:1<<max_level
		for p in 1:k
			# Check if our hier function is supported on this interval
			if (i-1)*h < left_point || i*h > right_point
				push!(ans, 0)
				continue
			end

			# Obtain the polynomial array for the position basis function
			unit_leg = shift_polyarray(
							scale_polyarray(leg_coeffs[p][1:k], 1//2), 1//2)
			pos_element = shift_polyarray(
								scale_polyarray(unit_leg, h), (i-1)*h)
			pos_element *= sqrt(1<<(max_level+1)) # (2.0)^((max_level+1)/2)

			# Perform symbolic integration depending on which
			# interval the position basis function overlaps with the
			# DG basis element
			if (i-1)*h < mid_point
				val = integrate_polyarray(vl, pos_element; a=(i-1)*h, b=i*h)
			else
				val = integrate_polyarray(vr, pos_element; a=(i-1)*h, b=i*h)
			end
			push!(ans, val)
		end
	end
	return ans
end


# Construct the change of basis matrix from hierarchical
# to position basis
function hier2pos(k::Int, max_level::Int; abs_tol=ABS_TOL)
	j = 1
	I = Int[]
	J = Int[]
	V = Float64[]
	for level in 0:max_level
		for cell in 1:(1<<max(0, level-1))
			for mode in 1:k
				ans = pos_vcoeffs_DG(k, max_level, v(k, level, cell, mode))
				for i in 1:length(ans)
					if abs(ans[i]) > abs_tol
						push!(I, i)
						push!(J, j)
						push!(V, ans[i])
					end
				end
				j += 1
			end
		end
	end
	return sparse(I, J, V, k * (1<<max_level), k * (1<<max_level), +)
end

function pos2hier(k::Int, max_level::Int; abs_tol=ABS_TOL)
	return hier2pos(k, max_level; abs_tol=abs_tol)'
end
