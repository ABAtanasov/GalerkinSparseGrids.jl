#------------------------------------------------------------
#
# Constructing the matrix elements of the derivative
# operator in the 1-D DG Basis
#
#------------------------------------------------------------

# Efficiency criticality: LOW
# Computations only performed once

# Accuracy criticality: HIGH
# Critical for accurate PDE evolution


#------------------------------------------------------
# 1-D Derivatives
#------------------------------------------------------

function dLegendreP(k,x)
    k<=K_max || throw(DomainError())
    return array2poly(symbolic_diff(leg_coeffs[k+1]),x)
end


function dh(k,f_number,x)
    f_number<=k || throw(DomainError())
    return array2poly(symbolic_diff((dg_coeffs[k])[f_number]),x)
end


#------------------------------------------------------
# Shifted and scaled derivatives: v
#------------------------------------------------------

function dv(k::Int, level::Int, place::Int, f_number::Int, x::Real)
	if level==0 # At the base level, we are not discontinuous, and we simply
                # use Legendre polynomials up to degree k-1 as a basis
        return 2*dLegendreP(f_number-1,2*x-1)*sqrt(2.0)
	else
        return (1<<(level))*dh(k, f_number, (1<<level)*x - (2*place-1)) * (2.0)^(level/2)
        # Otherwise we use an appropriately shifted, scaled, and normalized
		# DG function
	end
end

function dv(k::Int, level::Int, place::Int, f_number::Int)
    return (xs::Real -> dv(k,level,place,f_number,xs))
end


#------------------------------------------------------
# Numerical Inner Product
#------------------------------------------------------

# TODO: There is a way to do this all symbolically, with no use for numerics
function inner_product1D(f::Function, g::Function, lvl::Int, place::Int; 
								rel_tol = REL_TOL, abs_tol = ABS_TOL, max_evals=MAX_EVALS)
    xmin = (place-1)/(1<<(pos(lvl-1)))
	xmax = (place)/(1<<(pos(lvl-1)))
	h = (x-> f(x)*g(x))
    (val, err) = hquadrature(h, xmin, xmax; reltol=rel_tol, abstol=abs_tol, maxevals=max_evals)
	return val 
end

#------------------------------------------------------
# Taking Derivative of Array Representing Polynomial
#------------------------------------------------------

function symbolic_diff{T<:Real}(v::Array{T})
	n=length(v)
    k=div(n,2)
	ans = zeros(T,n)
	for i in 1:n
		if i<k
			ans[i] = i*v[i+1]
		elseif i > k && i<2k
			ans[i] = (i-k) * v[i+1]
		else
			ans[i]=0
		end
	end
	return ans
end

#------------------------------------------------------
# < f_1 | D | f_2 > matrix elements
#------------------------------------------------------

function legendreDlegendre(f_number1::Int, f_number2::Int)
	return inner_product(leg_coeffs[f_number1],symbolic_diff(leg_coeffs[f_number2]))
end 


function hDh(k::Int, f_number1::Int, f_number2::Int)
	return inner_product(dg_coeffs[k][f_number1], symbolic_diff(dg_coeffs[k][f_number2]))
end

function vDv(k::Int, lvl1::Int, place1::Int, f_number1::Int, lvl2::Int, place2::Int, f_number2::Int;
					rel_tol = REL_TOL, abs_tol = ABS_TOL, max_evals=MAX_EVALS)
	if lvl1 == lvl2
		if lvl1 == 0
			return 2*legendreDlegendre(f_number1, f_number2)
		end
		if place1 == place2
			return hDh(k, f_number1, f_number2)*(1<<pos(lvl1))
		end
		return 0.0
	end
	if lvl1 < lvl2
		if lvl1 == 0
            return inner_product1D(v(k,0,1,f_number1), dv(k,lvl2,place2,f_number2), lvl2, place2;
									rel_tol = rel_tol, abs_tol = abs_tol, max_evals=max_evals)
		end
		if (1<<(lvl2-lvl1))*place1 >= place2 && (1<<(lvl2-lvl1))*(place1-1) < place2
            return inner_product1D(v(k,lvl1,place1,f_number1), dv(k,lvl2,place2,f_number2), lvl2, place2;
									rel_tol = rel_tol, abs_tol = abs_tol, max_evals=max_evals)
		end
		return 0.0
    end
	return 0.0

end

#------------------------------------------------------
# < f_1 | D | f_2 > for a specific f_2 over all f_1
#------------------------------------------------------

function diff_basis_DG(k::Int, level::Int, place::Int, f_number::Int)
    dcoeffs = Array(Float64, (level+1,k))
    p = place
    for l in level:-1:0
        for f_n in 1:k
            dcoeffs[l+1,f_n]=vDv(k, l, p, f_n, level, place, f_number)
        end
        p = Int(ceil(p/2))
    end
    return dcoeffs
end

