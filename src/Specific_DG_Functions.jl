
const K_max = 10;

#Going from Arrays to Polynomials

@fastmath function array2poly{T<:Real}(v::Array{T},x::Real)
    if abs(x)>1 
		return zero(T)
	end
	n=length(v)
    k=div(n,2)
	s = zero(T)
	xi = one(T)
	@inbounds for i in k:-1:1
		# TODO: use Horner's method (Done)
		s *= x
		s += v[i] + flipsign(v[i+k],x)
	end
	return s
end

function array2poly{T<:Real}(v::Array{T})
    return (x-> array2poly(v,x))
end

#precompute Legendre

const Leg_coeffs=legendre(K_max)

function LegendreP(k,x)
    k<=K_max || throw(DomainError())
    return array2poly(Leg_coeffs[k+1],x)
end

function LegendreP(k)
    k<=K_max || throw(DomainError())
    return array2poly(Leg_coeffs[k+1])
end

#precompute DG functions

const DG_coeffs=Array(Array{Array{Float64,1},1}, K_max)
#TODO Make this a 2D array 
for i in 1:K_max
    DG_coeffs[i] = DG_Basis(i)
end

function h(k,f_number,x)
    f_number<=k || throw(DomainError())
    return array2poly((DG_coeffs[k])[f_number],x)
end

function h(k,f_number)
    f_number<=k || throw(DomainError())
    return array2poly((DG_coeffs[k])[f_number])
end

# Need to ask Erik some questions about precomputation
