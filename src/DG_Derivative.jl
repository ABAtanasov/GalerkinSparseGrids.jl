
#------------------------------------------------------
# 1-D Derivatives
#------------------------------------------------------

function dLegendreP(k,x)
    k<=K_max || throw(DomainError())
    return array2poly(symbolic_diff(Leg_coeffs[k+1]),x)
end

function dLegendreP(k)
    k<=K_max || throw(DomainError())
    return array2poly(symbolic_diff(Leg_coeffs[k+1]))
end

function dh(k,f_number,x)
    f_number<=k || throw(DomainError())
    return array2poly(symbolic_diff((DG_coeffs[k])[f_number]),x)
end

function dh(k,f_number)
    f_number<=k || throw(DomainError())
    return array2poly(symbolic_diff((DG_coeffs[k])[f_number]))
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
# Multidimensional derivatives
#------------------------------------------------------

function dV{D,T<:Real}(i::Int, k::Int, level::NTuple{D,Int}, 
    place::CartesianIndex{D}, f_number::CartesianIndex{D}, xs::AbstractArray{T})

    ans = one(eltype(xs))
    for j = 1:D
        if j == i
            ans *= dv(k, level[j], place[j], f_number[j], xs[j])
        else
            ans *= v(k, level[j], place[j], f_number[j], xs[j])
        end
    end
    return ans
end

function dV{D}(i, k, level::NTuple{D,Int}, 
                place::CartesianIndex{D}, f_number::CartesianIndex{D})
    return (xs::AbstractArray{Real} -> dV(i, k, level, place, f_number, xs))
end

#------------------------------------------------------
# Numerical Inner Product
#------------------------------------------------------

function inner_product1D(f::Function, g::Function, lvl::Int, place::Int)
    xmin = (place-1)/(1<<(pos(lvl-1)))
	xmax = (place)/(1<<(pos(lvl-1)))
	h = (x-> f(x)*g(x))
    (val, err) = hquadrature(h, xmin, xmax; reltol=1e-8, abstol=1e-8, maxevals=0)
	return val 
end
#There may be a way to do this all symbolically, with no use for numerics

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
#working

#------------------------------------------------------
# < f_1 | D | f_2 > matrix elements
#------------------------------------------------------

function legendreDlegendre(f_number1::Int, f_number2::Int)
	return inner_product(Leg_coeffs[f_number1],symbolic_diff(Leg_coeffs[f_number2]))
end 


function hDh(k::Int, f_number1::Int, f_number2::Int)
	return inner_product(DG_coeffs[k][f_number1], symbolic_diff(DG_coeffs[k][f_number2]))
end

function vDv(k::Int, lvl1::Int, place1::Int, f_number1::Int, lvl2::Int, place2::Int, f_number2::Int)
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
            return inner_product1D(v(k,0,1,f_number1), dv(k,lvl2,place2,f_number2), lvl2, place2)
		end
		if (1<<(lvl2-lvl1))*place1 >= place2 && (1<<(lvl2-lvl1))*(place1-1) < place2
            #@show (lvl1, place1, lvl2, place2)
            return inner_product1D(v(k,lvl1,place1,f_number1), dv(k,lvl2,place2,f_number2), lvl2, place2)
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


# The following code is not efficient for PDEsteps 
# and is replaced in the vMethods by appropriate sparse matrices

#
# #------------------------------------------------------
# # Modifying a coefficient set to get derivative coeffs
# # Asscociated with a specific basis function V
# #------------------------------------------------------
#
# function dchange{D,T<:Real}(k::Int,i::Int,dcoeffs::Dict{CartesianIndex{D}, Array{Array{Float64},D}},c::T,
#     lvl::NTuple{D,Int}, place::CartesianIndex{D}, f_number::CartesianIndex{D})
#     p = place[i]
#     for l in lvl[i]:-1:0
#         lvl1= ntuple(j-> j==i?l+1:(lvl[j]+1) ,D)
#         place1=ntuple(j-> j==i?p:place[j] ,D)
#         for f_n in 1:k
#             f_number1=ntuple(j-> j==i?f_n:f_number[j] ,D)
#             dcoeffs[CartesianIndex{D}(lvl1)][CartesianIndex{D}(place1)][CartesianIndex{D}(f_number1)] += c*precomputed_diffs[(k,lvl[i],place[i],f_number[i])][l+1,f_n]
#
#         end
#         p = Int(ceil(p/2))
#     end
# end
#
# #------------------------------------------------------
# # Finally, getting derivative coeffs from a coeff set
# #------------------------------------------------------
#
# function diff_coefficients_DG{D}(i::Int, k::Int,
#                         coeffs::Dict{CartesianIndex{D}, Array{Array{Float64},D}})
#
#     dcoeffs=deepcopy(coeffs)
#     f_numbers= ntuple(i-> k, D)
#
#     for key in keys(coeffs)
#         ks = ntuple(i -> 1<<pos(key[i]-2), D)
#         for place in CartesianRange(ks)
#             for f_number in CartesianRange(f_numbers)
#                 dcoeffs[key][place][f_number]=0
#             end
#         end
#     end
#
#     for key in keys(coeffs)
#         level = ntuple(i-> key[i]-1,D)
#         ks = ntuple(i -> 1<<pos(level[i]-1), D)
#         for place in CartesianRange(ks)
#             for f_number in CartesianRange(f_numbers)
#                 c=coeffs[key][place][f_number]
#                 dchange(k,i,dcoeffs,c,level,place,f_number)
#             end
#         end
#     end
#     return dcoeffs
# end
#
#
#
#
#
#
