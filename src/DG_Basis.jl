#------------------------------------------------------------
#
# Constructing the Hierarchical Discontinuous Galerkin Basis
#
#------------------------------------------------------------

# Efficiency criticality: LOW
# Computations only performed once


#------------------------------------------------------
# Defining the Inner Product
#------------------------------------------------------



# The coordinate convention here is to have vectors of length n even
# The first half has component i is the respective coefficient of x^i 
# THe second half has component n/2 + i is the coefficient of sgn(x)*x^i

function product_matrix(i::Int, j::Int, n::Int) 
    # This is the inner product of <x^i, x^j> when i,j < n/2 
    # or the appropriate reflection when they are >=
    # Exact expressions only. No NIntegrations
    k= Int(round(n/2))
    if i < k && j<k
        return (1 + (-1)^(i + j))/(1 + i + j) 
    elseif i >= k && j< k 
		return (1 - (-1)^((i-k) + (j)))/(1 + (i-k) + (j)) 
	elseif i < k && j >= k
		return product_matrix(j,i,n)
    else
		return (1 + (-1)^((i-k) + (j-k)))/(1 + (i-k) + (j-k)) 
    end
end

function inner_product{T<:Real}(v1::Array{T},v2::Array{T}) 
	value=zero(eltype(v1)) #Get 0 of the same type as v1
	n=length(v1)
	for i in 1:n
		for j in 1:n
            if v1[i]==0 || v2[j]==0
                continue 
            else
                value += product_matrix(i-1,j-1,n)*v1[i]*v2[j]
            end
		end
	end
    return value
end

function inner_product{T<:Real}(v1::Array{T},j::Int) #we consider x^(j-1)
	value=zero(eltype(v1)) #Get 0 of the same type as v1
	n=length(v1)
	for i in 1:n
        value += product_matrix(i-1,j-1,n)*v1[i] 
		#note we need to shift appropriately
	end
    return value
end

#------------------------------------------------------
# Defining Gram-Schmidt and forming Legendre Polynomials
#------------------------------------------------------

# The gram schmidt process on a set of vectors
# I think using an array of arrays is easiest here
function gram_schmidt{T<:Real}(Q_initial::Array{Array{T,1},1}) 
    n = length(Q_initial[1])
    k = Int(round(n/2)) 
	# n = 2k, k polys of degree up to k-1 and then k polys according to f
    Q_final = deepcopy(Q_initial) 
	# we'll return a modified version of the initial without modifying Q_initial
    for i = 1:k    
        for j = 1:i-1 
            proj = inner_product(Q_initial[i],Q_final[j])/inner_product(Q_final[j], Q_final[j]) 
			#This gives us the number corresponding to the projection along Q_final[j]
            Q_final[i] -= proj * Q_final[j]   
			#This subtracts it out from Q_final[i] to orthogonalize
        end
        Q_final[i] /= sqrt(inner_product(Q_final[i], Q_final[i])) 
		#Normalize every time
    end
    return Q_final
end

#We can now form the Legendre polynomials
function legendre(k::Int)
    Q = [[i==j?1.0:0.0 for i = 1:2*(k+1)] for j = 1:(k+1)]
	#start with just the x^i basis for i = 0 to k
    Q = gram_schmidt(Q)
	#Orthgonalize, and now we have the Legendre polynomials
    return Q
end
#working 

# I want to be able to do this completely analytically, using fractions
# so that I can avoid any numerical error in this process


#------------------------------------------------------
# Making a function have the first k moments vanish
#------------------------------------------------------

# # For just one function, v:
# function orthogonalize_1{T<:Real}(v::Array{T})     #modifier of v
#     n= length(v)
#     k= Int(round(n/2))
#     legendre_polys=legendre(k-1)
# 	#We need an orthogonal basis spanning the first k polys: x^0 to x^(k-1)
#     for j in 1:k
#         v[j]-= inner_product(v,legendre_polys[j])/product_matrix(j-1,j-1,n)
# 		#project out v by each element of this orthogonal basis
#     end
#     return v
# end
# #working

# For a basis of functions:
function orthogonalize_1{T<:Real}(Q_initial::Array{Array{T,1},1})
    n = length(Q_initial[1])
    k = Int(round(n/2)) 
	# n = 2k, because basis includes the f(j,x) functions in addition to x^j
    Q_final = deepcopy(Q_initial)
	# As before
    legendre_polys=legendre(k-1) 
	# We need an orthogonal basis for 1..x^(k-1) for the projection to work
    for i = 1:k
        for j in 1:k
            proj =  inner_product(Q_initial[i],legendre_polys[j])/
                        inner_product(legendre_polys[j],legendre_polys[j])
            Q_final[i] -= proj*legendre_polys[j] #subtract projection
			#project out v[i] by each element of this orthogonal basis for x^j
        end
    end
    return Q_final
end
#working

#------------------------------------------------------
# Make some functions have higher vanishing moments
#------------------------------------------------------

# This function will make k-1 of the basis vectors orth to x^k, 
# k-2 to x^k+1  all the way to 1 vector orth to x^2k-2
function orthogonalize_2{T<:Real}(Q_initial::Array{Array{T,1},1})
    n = length(Q_initial[1])
    k = Int(round(n/2))
    Q_final = deepcopy(Q_initial)
    for i = 1:k-1
        fi=copy(Q_initial[i]) 
		# We assume this isn't perp to x^(k+i-1) and subtract it from the next ones. 
		# We know that it isn't perp by parity considerations on the degree! :) 
		# (Thank goodness, cuz otherwise each time, we'd have to rearrange the basis
		# until we found a non-perp function and subtract by that one)
        for j = i+1:k
            a= inner_product(Q_final[j],k+i)/inner_product(fi, k+i)
            Q_final[j] -= a * fi
        end
    end 
    return Q_final
end


# Standard gram-schmidt process, starting at the end 
# (this is important, because the last function is orthogonal to a lot of higher 
# polynomials, and we don't want to do anything other than normalize it, 
# with similar reasoning for the penultimate, etc. functions) 
function gram_schmidt_rev{T<:Real}(Q_initial::Array{Array{T,1},1}) #I think an array of arrays is easiest here
    n = length(Q_initial[1])
    k = Int(round(n/2))
    Q_final = [[0.0 for i in 1:n] for j in 1:k]
    for i = k:-1:1    #note we're going in reverse 
        fi = copy(Q_initial[i])  #we won't modify the original, so we copy
        Q_final[i] = fi          #start with Q_final = f_i, and we'll subtract projections
        for j = i+1:k #because we're going in reverse
            proj = inner_product(fi,Q_final[j])/inner_product(Q_final[j], Q_final[j]) 
            #find the projection
            Q_final[i] -= proj * Q_final[j]    #project out direction j
        end
        Q_final[i] /= sqrt(inner_product(Q_final[i], Q_final[i])) #now normalize
    end
    return Q_final
end
#working


#------------------------------------------------------
# All together:
#------------------------------------------------------
function dg_basis(k::Int)
    Q = [[j==(i-k)?1.0:0.0 for i in 1:2*k] for j in 1:k]
    Q = orthogonalize_1(Q)
    Q = orthogonalize_2(Q) 
    Q = gram_schmidt_rev(Q)
    return Q
end