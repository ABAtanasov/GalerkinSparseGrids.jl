#------------------------------------------------------------
#
# Using boundary terms and integration (summation) by parts
# to construct more accurate derivative matrices for time
# evolution using Lax-Friedricks-type fluxes
#
# All derivatives are computed here in pos basis first
# Then conjugated into hier basis
#
#------------------------------------------------------------

# Efficiency criticality: LOW
# Computations only performed once

# Accuracy criticality: HIGH
# Critical for accurate PDE evolution

function leg(f_number::Int, x::Real)
    return sqrt(2.0)*LegendreP(f_number-1, 2*x-1) 
end

function basis(level::Int, place::Int, f_number::Int, x::Real)
    return leg(f_number, (1<<level)*x - (place-1)) * (2.0)^(level/2)
end

function basis(level::Int, place::Int, f_number::Int)
    return x->basis(level, place, f_number, x)
end

function dleg(f_number::Int, x::Real)
    return sqrt(2.0)*dLegendreP(f_number-1, 2*x-1) *2
end

function dbasis(level::Int, place::Int, f_number::Int, x::Real)
    return dleg(f_number, (1<<level)*x - (place-1)) * (2.0)^(level/2) * (1<<level)
end

function dbasis(level::Int, place::Int, f_number::Int)
    return x->dbasis(level, place, f_number, x)
end




function get_vcoeffs(k::Int, level::Int, f::Function; rel_tol = REL_TOL, abs_tol=ABS_TOL, max_evals=MAX_EVALS)
    vcoeffs = Array(Float64,(1<<level)*(k))
    i=1
    for place in 1:(1<<level)
        for f_number in 1:k
            vcoeffs[i] = hquadrature(x->(basis(level,place,f_number,x)*f(x)),
                            (place-1)/(1<<level), place/(1<<level), abstol=abs_tol)[1]
            i+=1
        end
    end
    return vcoeffs
end


# Not an efficient way to reconstrucy
# Dictionary lookup is better
function reconstruct_vcoeffs(k::Int, level::Int, vcoeffs::Array{Float64,1}, x::Real)
    value = 0.0
    i=1
    for place in 1:(1<<level)
        for f_number in 1:k
            value += vcoeffs[i]*basis(level,place,f_number,x)
            i+=1
        end
    end
    return value
end


function legvDv(level, place1, f_number1, place2, f_number2; rel_tol = REL_TOL, abs_tol=ABS_TOL, max_evals=MAX_EVALS)
    if place1 == place2
        return hquadrature(x->(basis(level,place1,f_number1,x)*dbasis(level,place2,f_number2,x)),
            (place1-1)/(1<<level), place1/(1<<level); abstol=abs_tol)[1]
    end
    return 0.0
end

function D_matrix(k::Int, level::Int)
    mat = spzeros((1<<level)*(k), (1<<level)*(k))
    i=1
    for place1 in 1:(1<<level)
        for f_number1 in 1:k
            j=1
            for place2 in 1:(1<<level)
                for f_number2 in 1:k
                    ans = legvDv(level, place1, f_number1, place2, f_number2)
                    if abs(ans) > 1.0e-15
						mat[i,j] = ans
					end
                    j+=1
                end
            end
            i+=1
        end
    end
    return mat
end



#------------------------------------------------------
# Lax-Friedrichs flux matrix element on Legendre basis
# This is currently supported only at alpha = 0
#------------------------------------------------------

# Periodic boundary:

function periodic_legvLFv(alpha::Real, level::Int, place1::Int, f_number1::Int, place2::Int, f_number2::Int)
    point1 = (place2-1)/(1<<level)
    point2 = (place2)/(1<<level)
    
    LF1 = 0.5* (basis(level,place1,f_number1,point1-5.0e-16)+basis(level,place1,f_number1,point1+5.0e-16)) + 
		   alpha *  (basis(level,place1,f_number1,point1+5.0e-16)-basis(level,place1,f_number1,point1-5.0e-16))
    LF2 = 0.5* (basis(level,place1,f_number1,point2-5.0e-16)+basis(level,place1,f_number1,point2+5.0e-16)) +
		   alpha *  (basis(level,place1,f_number1,point2+5.0e-16)-basis(level,place1,f_number1,point2-5.0e-16))
    
    if place2==(1<<level)
        LF2 = 0.5* (basis(level,place1,f_number1,point2-5.0e-16)+basis(level,place1,f_number1,0.0+5.0e-16)) +
		   alpha *  (basis(level,place1,f_number1,0.0+5.0e-16)-basis(level,place1,f_number1,point2-5.0e-16))
    end
    if place2==1
		LF1 = 0.5* (basis(level,place1,f_number1,1.0-5.0e-16)+basis(level,place1,f_number1,point1+5.0e-16)) + 
				   alpha *  (basis(level,place1,f_number1,point1+5.0e-16)-basis(level,place1,f_number1,1.0-5.0e-16))
    end

    val1 = basis(level,place2,f_number2,point1+5.0e-16)
    val2 = basis(level,place2,f_number2,point2-5.0e-16)
    
    return LF2 * val2 - LF1 * val1
end



function periodic_LF_matrix(alpha::Real, k::Int, level::Int)
    mat = spzeros((1<<level)*(k), (1<<level)*(k))
    i=1
    for place1 in 1:(1<<level)
        for f_number1 in 1:k
            j=1
            for place2 in 1:(1<<level)
                for f_number2 in 1:k
                    ans = periodic_legvLFv(alpha, level, place1, f_number1, place2, f_number2)
                    if abs(ans)<1.0e-15
                        j+=1
                        continue
                    end
                    mat[i,j]=ans
                    j+=1
                end
            end
            i+=1
        end
    end
    return mat
end

#------------------------------------------------------
# We can precompute the full matrix for this 
# boundary term
#------------------------------------------------------


function hier2pos(k::Int, max_level::Int; abs_tol = ABS_TOL)
    # mat = spzeros((1<<max_level)*(k), (1<<max_level)*(k))
    j=1
	I = Int[]
	J = Int[]
	K = Float64[]
	
    for level in 0:max_level
        for place in 1:(1<<pos(level-1))
            for f_number in 1:k
				# Issue is right here:
                ans = get_vcoeffs(k,max_level,x->v(k,level,place,f_number,x); abs_tol = abs_tol)
				
				for i in 1:length(ans)
                    if abs(ans[i]) > 1e-15
						push!(I, i)
						push!(J, j)
						push!(K, ans[i])
					# 	mat[i,j]=ans[i]
					# else
					# 	mat[i,j]=0
					end
                end
                j+=1
            end
        end
    end  
	mat = sparse(I, J, K, (1<<max_level)*(k), (1<<max_level)*(k), +)
    return mat
end


#hier2pos = Array(SparseMatrixCSC{Float64,Int64}, (KMAX, LMAX))


# # NOTE here this precomputation is done with ABS_TOL default.
# # Will want to change in the future
# for k in 1:KMAX
# 	for level in 1:LMAX
# 		hier2pos[k, level] = hier2pos_precompute(k, level)
# 	end
# end


# By default, use alpha = 0
function periodic_pos_DLF_Matrix(alpha::Real, k::Int, max_level::Int)
	A = -D_matrix(k,max_level) + periodic_LF_matrix(alpha, k,max_level)
	return A'
	
end


function periodic_hier_DLF_Matrix(alpha::Real, k::Int, max_level::Int)
	Q = hier2pos(k,max_level)
	A = periodic_pos_DLF_Matrix(alpha, k, max_level)
	return *(Q', *(A, Q))
	
end

