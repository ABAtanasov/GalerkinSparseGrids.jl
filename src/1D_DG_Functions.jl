#------------------------------------------------------------
#
# Evaluating the 1D Galerkin Basis functions
#
#------------------------------------------------------------

const K_max = 10;

# Efficiency criticality: HIGH
# Central to the efficiency of the package

# Going from arrays to polynomials

@fastmath function array2poly{T<:Real}(v::Array{T}, x::Real)
	if abs(x)>1 
		return zero(T)
	end
	n = length(v)
	k = div(n,2)
	s = zero(T)
	xi = one(T)
	@inbounds for i in k:-1:1
		# Using Horner's method 
		s *= x
		s += v[i] + flipsign(v[i+k], x)
	end
	return s
end

function array2poly{T<:Real}(v::Array{T})
	return (x-> array2poly(v,x))
end

#precompute Legendre

leg_coeffs=legendre(K_max)

function LegendreP(k, x)
	k<=K_max || throw(DomainError())
	return array2poly(leg_coeffs[k+1], x)
end

function LegendreP(k)
	k<=K_max || throw(DomainError())
	return array2poly(leg_coeffs[k+1])
end

#precomputing the DG functions

dg_coeffs = Array{Array{Array{Float64,1},1}}(K_max)
#TODO Make this a 2D array

for i in 1:K_max
	dg_coeffs[i] = dg_basis(i)
end


# This is the dg basis function corresponding to number mode
function h(k, mode, x)
	mode<=k || throw(DomainError())
	return array2poly((dg_coeffs[k])[mode],x)
end

function h(k, mode)
	mode<=k || throw(DomainError())
	return array2poly((dg_coeffs[k])[mode])
end
