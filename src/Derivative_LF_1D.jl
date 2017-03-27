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


function pos_vcoeffs_DG(k::Int, level::Int, f::Function;
						rel_tol = REL_TOL, abs_tol = ABS_TOL, max_evals = MAX_EVALS)
	vcoeffs = Array(Float64,(1<<level)*(k))
	i = 1
	for place in 1:(1<<level)
		for f_number in 1:k
			fcn			= x->(basis(level, place, f_number,x)*f(x))
			left_bound  = (place-1)/(1<<level)
			right_bound = (place)/(1<<level)
			vcoeffs[i]  = hquadrature(fcn, left_bound, right_bound; abstol=abs_tol)[1]
			i+=1
		end
	end
	return vcoeffs
end

function legvDv(level, place1, f_number1, place2, f_number2;
				rel_tol = REL_TOL, abs_tol=ABS_TOL, max_evals=MAX_EVALS)
	if place1 == place2
		fcn 		= x->(basis(level, place1,f_number1,x)*
							dbasis(level, place2, f_number2,x))
		left_bound  = (place1-1)/(1<<level)
		right_bound = (place1)/(1<<level)
		return hquadrature(fcn, left_bound, right_bound; abstol=abs_tol)[1]
	end
	return 0.0
end

function D_matrix(k::Int, level::Int)
	i = 1
	I = Int[]
	J = Int[]
	V = Float64[]
	for place1 in 1:(1<<level)
		for f_number1 in 1:k
			j = 1
			for place2 in 1:(1<<level)
				for f_number2 in 1:k
					ans = legvDv(level, place1, f_number1, place2, f_number2)
					if abs(ans) > 1.0e-15
						push!(I, i)
						push!(J, j)
						push!(V, ans)
					end
					j+=1
				end
			end
			i+=1
		end
	end
	return sparse(I, J, V, k * (1<<level), k * (1<<level), +)
end

#------------------------------------------------------
# Lax-Friedrichs flux matrix element on Legendre basis
# This is currently supported only at alpha = 0
#------------------------------------------------------

# Periodic boundary:

function periodic_legvLFv(level::Int, place1::Int, f_number1::Int,
									  place2::Int, f_number2::Int; alpha::Real = 0)
	point1 = (place2-1)/(1<<level)
	point2 = (place2)/(1<<level)
	tiny   = 5.0e-16

	left1  = basis(level, place1, f_number1, point1-tiny)
	right1 = basis(level, place1, f_number1, point1+tiny)
	left2  = basis(level, place1, f_number1, point2-tiny)
	right2 = basis(level, place1, f_number1, point2+tiny)

	(place2 == (1<<level)) && (right2 = basis(level, place1, f_number1, 0.0+tiny))
	(place2 == 1)		   && (left1  = basis(level, place1, f_number1, 1.0-tiny))

	LF1 = 0.5 *(left1 + right1) + alpha * (right1 - left1)
	LF2 = 0.5 * (left2 + right2) + alpha * (right2 - left2)

	val1 = basis(level, place2, f_number2, point1+tiny)
	val2 = basis(level, place2, f_number2, point2-tiny)

	# Net flux:
	return LF2 * val2 - LF1 * val1
end

function periodic_LF_matrix(k::Int, level::Int; alpha::Real = 0)
	i = 1
	I = Int[]
	J = Int[]
	V = Float64[]
	for place1 in 1:(1<<level)
		for f_number1 in 1:k
			j=1
			for place2 in 1:(1<<level)
				for f_number2 in 1:k
					ans = periodic_legvLFv(level, place1, f_number1,
												  place2, f_number2; alpha=alpha)
					if abs(ans)>1.0e-15
						push!(I, i)
						push!(J, j)
						push!(V, ans)
					end
					j+=1
				end
			end
			i+=1
		end
	end
	return sparse(I, J, V, k * (1<<level), k * (1<<level), +)
end

#------------------------------------------------------
# We can precompute the full matrix for this 
# boundary term
#------------------------------------------------------

function hier2pos(k::Int, max_level::Int; abs_tol = ABS_TOL)
	j = 1
	I = Int[]
	J = Int[]
	V = Float64[]
	for level in 0:max_level
		for place in 1:(1<<pos(level-1))
			for f_number in 1:k
				# Possible issue right here:
				ans = pos_vcoeffs_DG(k, max_level, x->v(k,level,place,f_number,x);
								  abs_tol = abs_tol)
				for i in 1:length(ans)
					if abs(ans[i]) > 1e-15
						push!(I, i)
						push!(J, j)
						push!(V, ans[i])
					end
				end
				j+=1
			end
		end
	end
	return sparse(I, J, V, k * (1<<max_level), k * (1<<max_level), +)
end


# By default, use alpha = 0
function periodic_pos_DLF_Matrix(k::Int, max_level::Int; alpha = 0)
	A = -D_matrix(k, max_level) + periodic_LF_matrix(k, max_level; alpha = alpha)
	return A'
end

function periodic_hier_DLF_Matrix(k::Int, max_level::Int; alpha = 0)
	Q = hier2pos(k, max_level)
	A = periodic_pos_DLF_Matrix(k, max_level; alpha = alpha)
	return *(Q', *(A, Q))
end

function periodic_DLF_Matrix(k::Int, max_level::Int; alpha::Real = 0, base = "hier")
	if base == "hier"
		return periodic_hier_DLF_Matrix(k, max_level; alpha = alpha)
	elseif base == "pos"
		return periodic_pos_DLF_Matrix(k, max_level; alpha = alpha)
	else
		throw(ArgumentError)
	end
end
