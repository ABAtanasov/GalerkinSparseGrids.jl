using ODE

function pos_wave_equation45(f0::Function, v0::Function, k::Int,level::Int, time0::Real, time1::Real)
	f0coeffs=get_vcoeffs(k,level, f0)
	v0coeffs=get_vcoeffs(k,level, v0)
	len = length(f0coeffs)
    laplac= *(periodic_pos_DLF_Matrix(0,k,level),periodic_pos_DLF_Matrix(0,k,level))
	RHS = spzeros(2*len, 2*len)
	for i in len+1:2*len
	    for j in 1:len
	        RHS[i,j] = laplac[i-len,j]
	        RHS[j,j+len] = 1.0
	    end
	end
    y0 = Array{Float64}([i<=len?f0coeffs[i]:v0coeffs[i-len] for i in 1:2*len])
	soln=ode45((t,x)->*(RHS,x), y0, [time0, time1])
	return soln
end

function hier_wave_equation45(f0::Function, v0::Function, k::Int,level::Int, time0::Real, time1::Real)
    f0coeffs=vhier_coefficients_DG(k,f0, (level+1,))
    v0coeffs=vhier_coefficients_DG(k,v0, (level+1,))
	len = length(f0coeffs)
    D_op = periodic_hier_DLF_Matrix(0,k,level)
    laplac= *(D_op,D_op)
	RHS = spzeros(2*len, 2*len)
	for i in len+1:2*len
	    for j in 1:len
	        RHS[i,j] = laplac[i-len,j]
	        RHS[j,j+len] = 1.0
	    end
	end
    y0 = Array{Float64}([i<=len?f0coeffs[i]:v0coeffs[i-len] for i in 1:2*len])
    soln=ode45((t,x)->*(RHS,x), y0, [time0,time1])
	return soln
end

function norm_squared{T<:Real}(coeffs::Array{T})
    sum = 0
    for i in coeffs
        sum+= i^2
    end
    return sum
end

function pos_energy_func(k, level, soln::Tuple{Array{Float64,1},Array{Array{Float64,1},1}})
    len = length(soln[1])
    num_coeffs = Int(round(length(soln[2][1])/2))
    
    times = copy(soln[1])
    energies = Array(Float64,len)
    
    D_op = periodic_pos_DLF_Matrix(0,k,level)

    for i in 1:len
        ux   = *(D_op,soln[2][i][1:num_coeffs])
        udot = soln[2][i][num_coeffs+1:end]
        energies[i] = norm_squared(ux)+norm_squared(udot)

    end
    return (times, energies)
end

function hier_energy_func(k, level, soln::Tuple{Array{Float64,1},Array{Array{Float64,1},1}})
    len = length(soln[1])
    num_coeffs = Int(round(length(soln[2][1])/2))
    
    times = copy(soln[1])
    energies = Array(Float64,len)
    
    D_op = periodic_hier_DLF_Matrix(0,k,level)

    for i in 1:len
        ux   = *(D_op,soln[2][i][1:num_coeffs])
        udot = soln[2][i][num_coeffs+1:end]
        energies[i] = norm_squared(ux)+norm_squared(udot)

    end
    return (times, energies)
end

# soln=wave_equation(x->sin(2*pi*x),x->0, 3,5,0,1)
#
# p = plotfunc1D(x->reconstruct_vcoeffs(3,5,soln[2][1], x[1]))
# anim = Animation()
#
# for i in 1:50:length(soln[2])
#                   p = plotfunc1D(x->reconstruct_vcoeffs(3,5,soln[2][i], x[1]))
# 				  xs=linspace(0,1,300)
# 				  ys=[1.0 for x in xs]
# 				  plot!(xs,ys)
# 				  ys=[-1.0 for x in xs]
# 				  plot!(xs,ys)
#                   frame(anim)
#               end
#
# gif(anim)