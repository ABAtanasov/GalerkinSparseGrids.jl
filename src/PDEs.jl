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

# Efficiency criticality: MED
# Most computations time is taken up by functions called from
# other script files

# Accuracy criticality: LOW
# The accuracy is limited only by functions in other scripts 
# and ODE.jl


function wave_evolve_1D(k::Int, max_level::Int,
							  f0::Function, v0::Function,
							  time0::Real, time1::Real; base = "hier", order = "45")
	if base == "pos"
		f0coeffs = pos_vcoeffs_DG(k, max_level, f0)
		v0coeffs = pos_vcoeffs_DG(k, max_level, v0)
	elseif base == "hier"
		f0coeffs = vcoeffs_DG(1, k, max_level, f0)
		v0coeffs = vcoeffs_DG(1, k, max_level, v0)
	else
		throw(ArgumentError(:base))
	end

	len = get_size(1, k, max_level)
	D_op = periodic_DLF_Matrix(k, max_level; base=base)
	laplac= *(D_op,D_op)
	RHS = spzeros(2*len, 2*len)

	for i in len+1:2*len
		for j in 1:len
			RHS[i,j] = laplac[i-len,j]
			RHS[j,j+len] = 1.0
		end
	end
	y0 = Array{Float64}([i<=len?f0coeffs[i]:v0coeffs[i-len] for i in 1:2*len])
	if order == "45"
		soln = ode45((t,x)->*(RHS,x), y0, [time0,time1])
	elseif order == "78"
		soln = ode78((t,x)->*(RHS,x), y0, [time0,time1])
	else
		throw(ArgumentError(:order))
	end
	return soln
end

function norm_squared{T<:Real}(coeffs::Array{T})
	sum = 0
	for i in coeffs
		sum += i^2
	end
	return sum
end

function energy_func_1D(k, level, 
						soln::Tuple{Array{Float64,1},Array{Array{Float64,1},1}}; 
						base = "hier")
	len			= length(soln[1])
	num_coeffs	= Int(round(length(soln[2][1])/2))
	times    	= copy(soln[1])
	energies	= Array(Float64, len)
	D_op		= periodic_DLF_Matrix(k, level; base=base)

	for i in 1:len
		ux   = *(D_op, soln[2][i][1:num_coeffs])
		udot = soln[2][i][num_coeffs+1:end]
		energies[i] = norm_squared(ux) + norm_squared(udot)

	end
	return (times, energies)
end

function wave_evolve(D::Int, k::Int, n::Int,
							  f0::Function, v0::Function,
							  time0::Real, time1::Real;
							  order = "45", scheme="sparse")
	println("0")
	f0coeffs = vcoeffs_DG(D, k, n, f0; scheme=scheme)
	v0coeffs = vcoeffs_DG(D, k, n, v0; scheme=scheme)
	println("1")
	laplac   = laplacian_matrix(D, k, n; scheme=scheme)
	println("2")
	len 	 = length(f0coeffs)
	I = Int[]
	J = Int[]
	V = Float64[]
	println("3")

	for i in len+1:2*len
		for j in 1:len
			push!(I, i)
			push!(J, j)
			push!(V, laplac[i - len, j])
		end
	end
	println("4")
	for j in 1:len
		push!(I, j)
		push!(J, j + len)
		push!(V, 1.0)
	end
	println("5")
	RHS = sparse(I, J, V, 2*len, 2*len, +)
	println("6")
	y0 = Array{Float64}([i<=len?f0coeffs[i]:v0coeffs[i-len] for i in 1:2*len])
	println("7")
	if order == "45"
		soln = ode45((t,x)->*(RHS,x), y0, [time0,time1])
	elseif order == "78"
		soln = ode78((t,x)->*(RHS,x), y0, [time0,time1])
	else
		throw(ArgumentError(:order))
	end
	println("8")
	return soln
end

function energy_func(D::Int, k::Int, n::Int, 
					 soln::Tuple{Array{Float64,1},Array{Array{Float64,1},1}};

	scheme = "sparse")
	len 		= length(soln[1])
	num_coeffs	= Int(round(length(soln[2][1])/2))
	times 		= copy(soln[1])
	energies	= Array(Float64, len)
	D_ops 		= grad_matrix(D, k, n; scheme=scheme)

	for i in 1:len
		ux   = [*(D_ops[j], soln[2][i][1:num_coeffs]) for j in 1:D]
		udot = soln[2][i][num_coeffs+1:end]
		energies[i] = sum([norm_squared(ux[j]) for j in 1:D])+norm_squared(udot)
	end
	return (times, energies)
end
