include("Specific_DG_Functions.jl")
include("DG_Derivative.jl")
include("DG_Methods.jl")

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



function get_coeffs(k::Int, level::Int, f::Function)
    coeffs = Array(Float64,(1<<level, k))
    for place in 1:(1<<level)
        for f_number in 1:k
            coeffs[place,f_number] = hquadrature(x->(basis(level,place,f_number,x)*f(x)),
            (place-1)/(1<<level), place/(1<<level), abstol=1.0e-10)[1]
        end
    end
    return coeffs
end

function get_vcoeffs(k::Int, level::Int, f::Function)
    vcoeffs = Array(Float64,(1<<level)*(k))
    i=1
    for place in 1:(1<<level)
        for f_number in 1:k
            vcoeffs[i] = hquadrature(x->(basis(level,place,f_number,x)*f(x)),
                            (place-1)/(1<<level), place/(1<<level), abstol=1.0e-10)[1]
            i+=1
        end
    end
    return vcoeffs
end



function reconstruct_coeffs(coeffs::Array{Float64,2}, x::Real)
    value = 0.0;
    level = Int(round(log2(length(coeffs[:,1]))))
    k = length(coeffs[1,:])
    for place in 1:(1<<level)
        for f_number in 1:k
            value += coeffs[place,f_number]*basis(level,place,f_number,x)
        end
    end
    return value
end

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


function legvDv(level, place1, f_number1, place2, f_number2)
    if place1 == place2
        return hquadrature(x->(basis(level,place1,f_number1,x)*dbasis(level,place2,f_number2,x)),
            (place1-1)/(1<<level), place1/(1<<level), abstol=1.0e-10)[1]
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
# Boundary jump (Delta) on the Legedre Basis against f
#------------------------------------------------------

function legfDeltav(level, f::Function, place2, f_number2)
    point1 = (place2-1)/(1<<level)
    point2 = (place2)/(1<<level)
    
    mean1 = 0.5* (f(point1-5.0e-16)+f(point1+5.0e-16))
    mean2 = 0.5* (f(point2-5.0e-16)+f(point2+5.0e-16))
    
    val1 = leg(f_number2, 0) * (2.0)^(level/2)
	#basis(level,place2,f_number2,point1)#+5.0e-16)
    val2 = leg(f_number2, 1) * (2.0)^(level/2)
	#basis(level,place2,f_number2,point2)#-5.0e-16)
    
    return mean2 * val2 - mean1 * val1
end


#------------------------------------------------------
# Lax-Friedrichs flux on the Legendre Basis against f
#------------------------------------------------------

function legfLFv(alpha, level, f::Function, place2, f_number2)
    point1 = (place2-1)/(1<<level)
    point2 = (place2)/(1<<level)
    
    LF1 = 0.5* (f(point1-5.0e-16)+f(point1+5.0e-16)) + 
				alpha * (1<<level) * (f(point1+5.0e-16)-f(point1-5.0e-16))
    LF2 = 0.5* (f(point2-5.0e-16)+f(point2+5.0e-16)) + 
				alpha * (1<<level) * (f(point1+5.0e-16)-f(point1-5.0e-16))
    
    val1 = basis(level,place2,f_number2,point1+5.0e-16)
    val2 = basis(level,place2,f_number2,point2-5.0e-16)
    
    return LF2 * val2 - LF1 * val1
end

#------------------------------------------------------
# Boundary jump (Delta) matrix elemtn on Legendre basis
#------------------------------------------------------

function legvDeltav(level::Int, place1::Int, f_number1::Int, place2::Int, f_number2::Int)
    point1 = (place2-1)/(1<<level)
    point2 = (place2)/(1<<level)
    
    mean1 = 0.5* (basis(level,place1,f_number1,point1-5.0e-16)+basis(level,place1,f_number1,point1+5.0e-16))
    mean2 = 0.5* (basis(level,place1,f_number1,point2-5.0e-16)+basis(level,place1,f_number1,point2+5.0e-16))
    
    val1 = leg(f_number2, 0) * (2.0)^(level/2)
    val2 = leg(f_number2, 1) * (2.0)^(level/2)
	
    if place2==(1<<level)
        mean2 = basis(level,place1,f_number1,point2-5.0e-16)
    end
    if place2==1
        mean1 = basis(level,place1,f_number1,point1+5.0e-16)
    end
    
    return mean2 * val2 - mean1 * val1
end

function legDeltavv(level, place1::Int, f_number1::Int, place2::Int, f_number2::Int)
    point1 = (place2-1)/(1<<level)
    point2 = (place2)/(1<<level)
    
    val1 = 1/2 * leg(f_number2, 0) * (2.0)^(level/2)
    val2 = 1/2 * leg(f_number2, 1) * (2.0)^(level/2)
    
	diff1 = (basis(level,place1,f_number1,point1+5.0e-16)-basis(level,place1,f_number1,point1-5.0e-16))
	diff2 = -(basis(level,place1,f_number1,point2+5.0e-16)-basis(level,place1,f_number1,point2-5.0e-16))
	
    # if place2==(1<<level)
#         diff2 = 0 #basis(level,place1,f_number1,point2-5.0e-16)
#     end
#     if place2==1
#         diff1 = 0 #basis(level,place1,f_number1,point1+5.0e-16)
#     end
    
    
    return  diff2 * val2 - diff1 * val1
end



#------------------------------------------------------
# Lax-Friedrichs flux matrix element on Legendre basis
#------------------------------------------------------

function legvLFv(alpha::Real, level::Int, place1::Int, f_number1::Int, place2::Int, f_number2::Int)
    point1 = (place2-1)/(1<<level)
    point2 = (place2)/(1<<level)
    
    LF1 = 0.5* (basis(level,place1,f_number1,point1-5.0e-16)+basis(level,place1,f_number1,point1+5.0e-16)) + 
		   alpha *  (basis(level,place1,f_number1,point1+5.0e-16)-basis(level,place1,f_number1,point1-5.0e-16))
    LF2 = 0.5* (basis(level,place1,f_number1,point2-5.0e-16)+basis(level,place1,f_number1,point2+5.0e-16)) +
		   alpha *  (basis(level,place1,f_number1,point2+5.0e-16)-basis(level,place1,f_number1,point2-5.0e-16))
    
    if place2==(1<<level)
        LF2 = basis(level,place1,f_number1,point2-5.0e-16)
    end
    if place2==1
        LF1 = basis(level,place1,f_number1,point1+5.0e-16)
    end

    val1 = basis(level,place2,f_number2,point1+5.0e-16)
    val2 = basis(level,place2,f_number2,point2-5.0e-16)
    
    return LF2 * val2 - LF1 * val1
end

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

#------------------------------------------------------
# We can precompute the full matrix for this 
# boundary term
#------------------------------------------------------

function Delta_matrix(k::Int, level::Int)
    mat = spzeros((1<<level)*(k), (1<<level)*(k))
    i=1
    for place1 in 1:(1<<level)
        for f_number1 in 1:k
            j=1
            for place2 in 1:(1<<level)
                for f_number2 in 1:k
                    ans = legvDeltav(level, place1, f_number1, place2, f_number2)
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

function Delta_matrix2(k::Int, level::Int)
    mat = spzeros((1<<level)*(k), (1<<level)*(k))
    i=1
    for place1 in 1:(1<<level)
        for f_number1 in 1:k
            j=1
            for place2 in 1:(1<<level)
                for f_number2 in 1:k
                    ans = legDeltavv(level, place1, f_number1, place2, f_number2)
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


function LF_matrix(alpha::Real, k::Int, level::Int)
    mat = spzeros((1<<level)*(k), (1<<level)*(k))
    i=1
    for place1 in 1:(1<<level)
        for f_number1 in 1:k
            j=1
            for place2 in 1:(1<<level)
                for f_number2 in 1:k
                    ans = legvLFv(alpha, level, place1, f_number1, place2, f_number2)
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


function hier2pos(k::Int, max_level::Int)
    mat = spzeros((1<<max_level)*(k), (1<<max_level)*(k))
    j=1
    for level in 0:max_level
        for place in 1:(1<<pos(level-1))
            for f_number in 1:k
                ans = get_vcoeffs(k,max_level,x->v(k,level,place,f_number,x))
                for i in 1:length(ans)
                    if abs(ans[i])<1.0e-15
                        continue
                    end
                    mat[i,j]=ans[i]
                end
                j+=1
            end
        end
    end  
    return mat
end

function pos_Derivative_Matrix(k::Int, max_level::Int)
	A = -D_matrix(k,max_level) + Delta_matrix(k,max_level)
	return A'
	
end

function pos_DLF_Matrix(alpha::Real, k::Int, max_level::Int)
	A = -D_matrix(k,max_level) + LF_matrix(alpha, k,max_level)
	return A'
	
end

function periodic_pos_DLF_Matrix(alpha::Real, k::Int, max_level::Int)
	A = -D_matrix(k,max_level) + periodic_LF_matrix(alpha, k,max_level)
	return A'
	
end


function hier_Derivative_Matrix(k::Int, max_level::Int)
	Q = hier2pos(k,max_level)
	A = -D_matrix(k,max_level) + Delta_matrix(k,max_level)
	return *(Q', *(A', Q))
	
end

function periodic_hier_DLF_Matrix(alpha::Real, k::Int, max_level::Int)
	Q = hier2pos(k,max_level)
	A = periodic_pos_DLF_Matrix(alpha, k, max_level)
	return *(Q', *(A, Q))
	
end
