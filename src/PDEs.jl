# -----------------------------------------------------------
#
# Using boundary terms and integration (summation) by parts
# to construct more accurate derivative matrices for time
# evolution using Lax-Friedricks-type fluxes
#
# All derivatives are computed here in pos basis first
# Then conjugated into hier basis
#
# -----------------------------------------------------------

# Efficiency criticality: MED
# Most computations time is taken up by functions called from
# other script files

# Accuracy criticality: LOW
# The accuracy is limited only by functions in other scripts
# and ODE.jl

# Construct the time evolution matrix given a laplacian and
# initial data given f0, v0
function wave_data(laplac::A, f0coeffs::Array{T, 1},
	v0coeffs::Array{T, 1}) where {A <: AbstractArray, T <: Real}

	len = length(f0coeffs)
	rows = rowvals(laplac)
	vals = nonzeros(laplac)

	I = Int[]
	J = Int[]
	V = T[]
	for col = 1:len
		for i in nzrange(laplac, col)
			row = rows[i]
			val = vals[i]
			push!(I,row+len)
			push!(J,col)
			push!(V,val)
		end
	end
	for i = 1:len
		push!(I, i)
		push!(J, i+len)
		push!(V, 1.0)
	end
	RHS = sparse(I, J, V, 2*len, 2*len, +)
	y0 = Array{T}(undef, [i<=len ? f0coeffs[i] : v0coeffs[i-len] for i in 1:2*len])
	return RHS, y0
end

# The reason we have a wave evolution in 1D is to test
# that the standard "position" and multiresolution "heirarchical"
# bases agree for wave PDE evolution
function wave_evolve_1D(k::Int, max_level::Int,
							f0::Function, v0::Function,
							time0::Real, time1::Real; basis="hier", order = "45",
							kwargs...)
	if basis == "pos"
		f0coeffs = pos_vcoeffs_DG(k, max_level, f0)
		v0coeffs = pos_vcoeffs_DG(k, max_level, v0)
	elseif basis == "hier"
		f0coeffs = vcoeffs_DG(1, k, max_level, f0)
		v0coeffs = vcoeffs_DG(1, k, max_level, v0)
	elseif basis == "nodal"
		throw(MethodError("Nodal basis wave evolution not yet implemented"))
	elseif basis == "point"
		m2p = modal2points_1D(k, max_level)
		f0coeffs = m2p*vcoeffs_DG(1, k, max_level, f0)
		v0coeffs = m2p*vcoeffs_DG(1, k, max_level, v0)
	else
		throw(ArgumentError(:basis))
	end

	D_op = periodic_DLF_matrix(k, max_level; basis=basis)
	laplac= *(D_op,D_op)

	RHS, y0 = wave_data(laplac, f0coeffs, v0coeffs)
	if order == "45"
		soln = ode45((t,x)->*(RHS,x), y0, [time0,time1]; kwargs...)
	elseif order == "78"
		soln = ode78((t,x)->*(RHS,x), y0, [time0,time1]; kwargs...)
	else
		throw(ArgumentError(:order))
	end
	return soln
end

function norm_squared(coeffs::Array{T},; basis = "hier") where T <: Real

	if basis=="hier" || basis=="pos"
		return sum(i^2 for i in coeffs)
	else
		throw(MethodError("Nodal and point basis norms not yet implemented"))
	end
end

function energy_func_1D(k, level, soln::Tuple{Array{T, 1}, Array{Array{T, 1}, 1}};
	basis = "hier") where T <: Real

	len			= length(soln[1])
	num_coeffs	= Int(round(length(soln[2][1])/2))
	times		= copy(soln[1])
	energies	= Array{T}(undef, len)
	D_op		= periodic_DLF_matrix(k, level; basis=basis)

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
							  order="45", scheme="sparse", kwargs...)

	f0coeffs = vcoeffs_DG(D, k, n, f0; scheme=scheme)
	v0coeffs = vcoeffs_DG(D, k, n, v0; scheme=scheme)
	laplac   = laplacian_matrix(D, k, n; scheme=scheme)

	RHS, y0 = wave_data(laplac, f0coeffs, v0coeffs)
	if order == "45"
		soln = ode45((t,x)->*(RHS,x), y0, [time0,time1]; kwargs...)
	elseif order == "78"
		soln = ode78((t,x)->*(RHS,x), y0, [time0,time1]; kwargs...)
	else
		throw(ArgumentError(:order))
	end
	return soln
end

function energy_func(D::Int, k::Int, n::Int, soln::Tuple{Array{T, 1}, Array{Array{T, 1}, 1}};
		scheme = "sparse") where T <: Real

	len 		= length(soln[1])
	num_coeffs	= Int(round(length(soln[2][1])/2))
	times 		= copy(soln[1])
	energies	= Array{T}(undef, len)
	D_ops		= grad_matrix(D, k, n; scheme=scheme)

	for i in 1:len
		ux   = [*(D_ops[j], soln[2][i][1:num_coeffs]) for j in 1:D]
		udot = soln[2][i][num_coeffs+1:end]
		energies[i] = sum([norm_squared(ux[j]) for j in 1:D])+norm_squared(udot)
	end
	return (times, energies)
end
