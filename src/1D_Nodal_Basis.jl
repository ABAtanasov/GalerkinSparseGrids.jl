# -----------------------------------------------------------
#
# Implementing a nodal basis for DG sparse grids 
#
# -----------------------------------------------------------

#= 

We begin by defining analogues of our DG functions
to satisfy the following conditions:

	1) There are a set of collocation points respecting the hierarchical structure
	2) The nodal basis functions at a given level are those that vanish
		on all nodal points at the above levels and vanish on all but one
		of the cells on the level below
	3) Third constraint?

To do this, we _restrict_ to the continuous part of the DG basis, and focus
on piecewise polynomials that vanish on their boundaries. We then pick the 
collocation points to be exactly at these boundaries.

For this, it is useful to define lagrange polynomials, and moreover 
use k=2, 3, 5. 
The k=2 case has our basis is exactly the hat function basis constructed 
previously.

For k=2, the coarsest level has nodal points at {0, 1}, 
	the next level (the first hat) has nodal point at 0.5
	the next level has two cells with points at 0.25, 0.75 etc.

For k=3 the coarsest level has nodal points at {0, 0.5, 1},
	the next level has two nodal functions with points at 0.25, 0.75,
	the next level has two cells and two nodes per cell at all odd multiples of 0.125 etc.

For k=5 the coarsest level has nodal points at {0, 0.25, 0.5, 0.75, 1}
	the next level has four nodal functions with points at all multiples of 0.125
	the next level has two cells and four nodes per cell at all odd multiples of 0.0625 etc.

Notice that k=2,3,5 work because at each cell of the multiresolution basis the space of
	picewise continuous polynomials is k-1 = 2^i for i=0,1,2 in this case.
	A construction for other k would be harder, if not impossible.

=#

# Lagrange polynomials interpolating at the coarsest level
function lag_nodal(k::Int, mode::Int, x::Real)
	(k < 2 || k > 5) && throw(ArgumentError("k is not between 2 and 5"))
	(mode > k || mode < 1) && throw(ArgumentError("mode is not in [1, k]"))
	if k == 2
		mode == 1 && return max(zero(x), x)
		mode == 2 && return max(zero(x), 1-x)
	elseif k == 3
		mode == 1 && return (x < 0 || x > 1) ? zero(x) : 2 * (x-1) * (x - 1/2)
		mode == 2 && return (x < 0 || x > 1) ? zero(x) : 2 * x * (x - 1/2)
		mode == 3 && return (x < 0 || x > 1) ? zero(x) : -4 * x * (x-1)
	elseif k == 4 
		throw(MethodError("k = 4 not implemented"))
	elseif k == 5
		mode == 1 && return (x > 1 || x < 0) ? zero(x) : 32/3 * (x - 1/4) * (x - 1/2) * (x - 3/4) * (x - 1)
		mode == 2 && return (x > 1 || x < 0) ? zero(x) : -128/3 * x * (x - 1) * (x - 1/2) * (x - 3/4)
		mode == 3 && return (x > 1 || x < 0) ? zero(x) : 64 * x * (x - 1) * (x - 1/4) * (x - 3/4)
		mode == 4 && return (x > 1 || x < 0) ? zero(x) : -128/3 * x * (x - 1) * (x - 1/2) * (x - 1/4)
		mode == 5 && return (x > 1 || x < 0) ? zero(x) : 32/3 * x * (x - 1/4) * (x - 1/2) * (x - 3/4)
	end
end

# Multiresolution functions for the levels below
function h_nodal(k::Int, mode::Int, x::Real)
	(k < 2 || k > 5) && throw(ArgumentError("k is not between 2 and 5"))
	(mode >= k) && return zero(x)
	#(abs(x) >= 1) && return zero(T)
	if k == 2
		return max(zero(x), 1-abs(x))
	elseif k == 3
		mode == 1 && return max(zero(x), -4 * x * (x + 1))
		mode == 2 && return max(zero(x), -4 * x * (x - 1))
	elseif k == 4
		throw(MethodError("k = 4 not implemented"))
	elseif k == 5
		mode == 1 && return (x < -1 || x > 0) ? zero(x) : -128/3 * x * (x + 1) * (x + 1/2) * (x + 1/4)
		mode == 2 && return (x < -1 || x > 0) ? zero(x) : -128/3 * x * (x + 1) * (x + 1/2) * (x + 3/4)
		mode == 3 && return (x > 1  || x < 0) ? zero(x) : -128/3 * x * (x - 1) * (x - 1/2) * (x - 3/4)
		mode == 4 && return (x > 1  || x < 0) ? zero(x) : -128/3 * x * (x - 1) * (x - 1/2) * (x - 1/4)
	end
end

# This yields our nodal basis:
function v_nodal(k::Int, level::Int, cell::Int, mode::Int, x::Real)
	if level==0
		return lag_nodal(k, mode, x)
	else
		return h_nodal(k, mode, (1<<level)*x - (2*cell-1))
	end
end

function v_nodal(k::Int, level::Int, cell::Int, mode::Int)
	return x->v_nodal(k, level, cell, mode, x)
end

# -----------------------------------------------------------
# Next we define the sparse matrices transforming between 
# the heirarchical, nodal, and point bases
# -----------------------------------------------------------

# Evaluation of a given nodal function along all relevant gridpoints
function eval_points_1D(k::Int, max_level::Int, level::Int, cell::Int, mode::Int)
	i = 1
	I = Int[]
	V = Float64[]

	# level 0:
	for x in 0:(1/(k-1)):1
		val = v_nodal(k, level, cell, mode, x)
		(val != zero(x)) && (push!(I, i); push!(V, val))
		i+=1
	end
	# subsequent levels:
	for l in 1:max_level
		for c in 1:1<<max(0, l-1)
			for node in 1/(2*(k-1)):1/(k-1):(1-1/(2*(k-1)))
				x = (c - 1 + node)/(1<<(l-1))
				val = v_nodal(k, level, cell, mode, x)
				(val != zero(x)) && (push!(I, i); push!(V, val))
				i += 1
			end
		end
	end
	return dropzeros!(sparsevec(I, V))
end

# Sparse matrix converting nodal basis element to its value
# on the collocation points
function nodal2points_1D(k::Int, max_level::Int)
	I = Int[]
	J = Int[]
	V = Float64[]
	j = 1
	for mode in 1:k
		coeffs = eval_points_1D(k, max_level, 0, 1, mode)
		for (i, val) in zip(findnz(coeffs)...)
			push!(I, i);
			push!(J, j)
			push!(V, val)
		end
		j += 1
	end
	for level in 1:max_level
		for cell in 1:1<<max(0, level-1)
			for mode in 1:(k-1)
				coeffs = eval_points_1D(k, max_level, level, cell, mode)
				for (i, val) in zip(findnz(coeffs)...)
					push!(I, i)
					push!(J, j)
					push!(V, val)
				end
				j += 1
			end
		end
	end
    return dropzeros!(sparse(I, J, V))
end

# Sparse matrix inverting the above construction
# Converting a set of values on the collocation points to the DG space
# in terms of the (for now continuous) nodal basis
# TODO: Allow for this to handle discontinuity across cells
function points2nodal_1D(k::Int, max_level::Int)
    Q = nodal2points_1D(k, max_level)
    return threshold(inv(full(Q)), abs_tol=1e-15)
end

# Converts a given nodal basis element to a sum in the standard 
# hieararchical DG scheme
function nodal2heir_1D(k::Int, level::Int, cell::Int, mode::Int; rel_tol=1e-10, abs_tol=1e-12, max_evals=1500)
	I = Int[]
	V = Float64[]
	i = 1
	
	for hier_level in 0:level
		hier_cells = 1<<max(0, hier_level-1)
		for hier_cell in 1:hier_cells
			for hier_mode in 1:k
				x_min = (hier_cell-1)/hier_cells
				x_max = (hier_cell)/hier_cells
				val = hquadrature(x->v_nodal(k, level, cell, mode, x)*
									v(k, hier_level, hier_cell, hier_mode, x),
									x_min, x_max,
									reltol=rel_tol, abstol=abs_tol, maxevals=max_evals)[1]
				(abs(val) > eps(Float64)) && (push!(I, i); push!(V, val))
				i += 1
			end
		end
	end
	return dropzeros!(sparsevec(I, V))
end

function nodal2heir_1D(k::Int, max_level::Int; rel_tol=1e-10, abs_tol=1e-12, max_evals=1500)
	I = Int[]
	J = Int[]
	V = Float64[]
	j = 1
	
	for nodal_level in 0:max_level
		nodal_level == 0 ? num_nodes = k : num_nodes = k-1
		nodal_cells = 1<<max(0, nodal_level-1)
		for nodal_cell in 1:nodal_cells
			for nodal_mode in 1:num_nodes
				coeffs = nodal2heir_1D(k, nodal_level, nodal_cell, nodal_mode;
										rel_tol=rel_tol,
										abs_tol=abs_tol,
										max_evals=max_evals)
				for (i, val) in zip(findnz(coeffs)...)
					push!(I, i)
					push!(J, j)
					push!(V, val)
				end
				j += 1
			end
		end
	end
	return sparse(I, J, V)
end

# Accounting for ambiguity in the definition of the hierarchical
# basis across regions
function eval_v(k, level, cell, mode, x)
	xl = max(x - 1e-16, 0)
	xr = min(x + 1e-16, 1)
	return 0.5 * (v(k, level, cell, mode, xl) + v(k, level, cell, mode, xr))
end

# Evaluates a given hierarchical basis function at
# all collocation points
function hier2points_1D(k::Int, max_level::Int, level::Int, cell::Int, mode::Int)
	i = 1
	I = Int[]
	V = Float64[]

	# level 0:
	for x in 0:(1/(k-1)):1
		val = eval_v(k, level, cell, mode, x)
		(abs(val) != zero(x)) && (push!(I, i); push!(V, val))
		i+=1
	end
	
	# subsequent levels:
	for l in 1:max_level
		for c in 1:1<<max(0, l-1)
			for node in 1/(2*(k-1)):1/(k-1):(1-1/(2*(k-1)))
				x = (c - 1 + node)/(1<<(l-1))
				val = eval_v(k, level, cell, mode, x)
				(abs(val) != 0) && (push!(I, i); push!(V, val))
				i += 1
			end
		end
	end
	return dropzeros!(sparsevec(I, V))
end

# Constructs the sparse matrix to convert from the 
# hierarchical basis to the point basis
function hier2points_1D(k::Int, max_level::Int)
	I = Int[]
	J = Int[]
	V = Float64[]
	j = 1
	for mode in 1:k
		coeffs = hier2points_1D(k, max_level, 0, 1, mode)
		for (i, val) in zip(findnz(coeffs)...)
			push!(I, i)
			push!(J, j)
			push!(V, val)
		end
		j += 1
	end
	for level in 0:max_level
		for cell in 1:1<<max(0, level-1)
			for mode in 1:k
				coeffs = hier2points_1D(k, max_level, level, cell, mode)
				for (i, val) in zip(findnz(coeffs)...)
					push!(I, i)
					push!(J, j)
					push!(V, val)
				end
				j += 1
			end
		end
	end
	return dropzeros!(sparse(I, J, V))
end

# Inverts the above construction by factoring through
# the nodal space
function points2hier_1D(k::Int, max_level::Int)
	p2n = points2nodal_1D(k, max_level)
	n2h = nodal2heir_1D(k, max_level)
	return threshold(n2h*p2n)
end



