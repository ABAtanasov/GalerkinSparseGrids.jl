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

# Efficiency criticality: MEDIUM
# Computations performed once, 
# but take very long in multidimensional setting

# Accuracy criticality: HIGH
# Critical for accurate PDE evolution

# After the integration by parts, this corresponds to 
# the remaining integral, ignoring the boundary term
function D_matrix(k::Int, level::Int)
	i = 1
	I = Int[]
	J = Int[]
	V = Float64[]
	for cell1 in 1:(1<<level)
		for mode1 in 1:k
			j = 1
			for cell2 in 1:(1<<level)
				for mode2 in 1:k
					ans = legvDv(level, cell1, mode1, cell2, mode2)
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
# For now, use alpha = 0 only
#------------------------------------------------------

# Periodic boundary term:
function periodic_legvLFv(level::Int, cell1::Int, mode1::Int,
									  cell2::Int, mode2::Int; alpha::Real = 0)
	point1 = (cell2-1)/(1<<level)
	point2 = (cell2)/(1<<level)
	tiny   = 5.0e-16

	left1  = basis(level, cell1, mode1, point1-tiny)
	right1 = basis(level, cell1, mode1, point1+tiny)
	left2  = basis(level, cell1, mode1, point2-tiny)
	right2 = basis(level, cell1, mode1, point2+tiny)

	(cell2 == (1<<level)) && (right2 = basis(level, cell1, mode1, 0.0+tiny))
	(cell2 == 1)		   && (left1  = basis(level, cell1, mode1, 1.0-tiny))

	LF1 = 0.5 *(left1 + right1) + alpha * (right1 - left1)
	LF2 = 0.5 * (left2 + right2) + alpha * (right2 - left2)

	val1 = basis(level, cell2, mode2, point1+tiny)
	val2 = basis(level, cell2, mode2, point2-tiny)

	# Net flux:
	return LF2 * val2 - LF1 * val1
end

function periodic_LF_matrix(k::Int, level::Int; alpha::Real = 0)
	i = 1
	I = Int[]
	J = Int[]
	V = Float64[]
	for cell1 in 1:(1<<level)
		for mode1 in 1:k
			j=1
			for cell2 in 1:(1<<level)
				for mode2 in 1:k
					ans = periodic_legvLFv(level, cell1, mode1,
												  cell2, mode2; alpha=alpha)
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

#----------------------------------------------------------
# We can now build the full derivative matrix  
# in both the position and hierarchical basis
# with some help from the hier2pos change of basis matrix
#----------------------------------------------------------

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
