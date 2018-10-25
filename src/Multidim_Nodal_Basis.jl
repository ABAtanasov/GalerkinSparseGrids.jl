# -----------------------------------------------------------
#
# Implementing a nodal basis for DG sparse grids,
# Multidimensional case
#
# -----------------------------------------------------------

#=

	We continue the work started in 1D_Nodal_Basis.jl by applying the
	tensor product construction together with the sparse cutoff
	to constract a way to translate between the modal basis used in
	linear wave evolution and the nodal/point basis necessary for
	efficiently handling nonlinearities in the time evolution equations

=#

# Check whether (level1, cell1) overlaps with (level2, cell2)
function relevant_cell_1D(k::Int, level1::Int, cell1::Int, level2::Int, cell2::Int)
    left1  = (cell1-1)/(1<<max(0, level1-1))
    right1 = cell1/(1<<max(0, level1-1))
    left2  = (cell2-1)/(1<<max(0, level2-1))
    right2 = cell2/(1<<max(0, level2-1))
    return (left1 <= left2 && right1 >= right2) || (left1 >= left2 && right1 <= right2)
end

# Multidimensional version of the above
function relevant_cell(k::Int, level1::CartesianIndex{D},
	 	cell1::CartesianIndex{D}, level2::CartesianIndex{D},
	 	cell2::CartesianIndex{D}) where D

	bool = true
	for d in 1:D
		bool &= relevant_cell_1D(k, level1[d]-1, cell1[d], level2[d]-1, cell2[d])
	end
	return bool
end

# Takes an input (level, cell, mode) and yielding the
# corresponding index in the vector of coefficients for the nodal/point basis
function get_point_index(k::Int, l::Int, c::Int, m::Int)
	if l == 0
		return m
	end
	return k + (k-1) * (1<<(l-1) - 1) + (k-1) * (c-1) + m
end

function get_modal_index(k::Int, l::Int, c::Int, m::Int)
	if l == 0
		return m
	end
	return k * (1<<(l-1) + (c-1) ) + m
end

function inner_loop_to_points(c2::Int, k::Int, l::CartesianIndex{D},
			level::CartesianIndex{D}, cell::CartesianIndex{D},
			mode::CartesianIndex{D}, I::Array{Int, 1}, V::Array{T, 1},
			base_indices::Array{Int, 1}, point_mat_1D::Array{T, 2}) where {D, T <: Real}

	cells::NTuple{D, Int} = ntuple(i -> 1<<max(0, l[i]-2), D)
	modes::NTuple{D, Int} = ntuple(i->l[i]==1 ? k : k-1, D)
	cellrange = CartesianIndex(cells)
	moderange = CartesianIndex(modes)
	for c in cellrange
		!relevant_cell(k, l, c, level, cell) && (c2 += prod(modes); continue)

		for m in moderange
			val = one(T)
			for d in 1:D
				index = get_point_index(k, l[d]-1, c[d], m[d])
				val *= point_mat_1D[index, base_indices[d]]
			end
			(abs(val) < eps(T)) && (c2 += 1; continue)
			push!(I, c2)
			push!(V, val)
			c2 += 1
		end
	end
	return c2
end

function inner_loop_to_modal(c2::Int, k::Int, l::CartesianIndex{D}, level::CartesianIndex{D},
			cell::CartesianIndex{D}, mode::CartesianIndex{D}, I::Array{Int, 1}, V::Array{T, 1},
			base_indices::Array{Int, 1}, point_mat_1D::Array{T, 2}) where {D, T <: Real}

	cells::NTuple{D, Int} = ntuple(i -> 1<<max(0, l[i]-2), D)
	modes::NTuple{D, Int} = ntuple(i -> k, D)
	cellrange = CartesianIndex(cells)
	moderange = CartesianIndex(modes)
	for c in cellrange
		!relevant_cell(k, l, c, level, cell) && (c2 += prod(modes); continue)

		for m in moderange
			val = one(T)
			for d in 1:D
				index = get_modal_index(k, l[d]-1, c[d], m[d])
				val *= point_mat_1D[index, base_indices[d]]
			end
			(abs(val) < eps(T)) && (c2 += 1; continue)
			push!(I, c2)
			push!(V, val)
			c2 += 1
		end
	end
	return c2
end



# Evaluates a given multidimensional basis element in the modal (DG Sparse) basis
# on all the points of the corresponding sparse grid
function eval_points(k::Int, n::Int, level::CartesianIndex{D},
				cell::CartesianIndex{D}, mode::CartesianIndex{D},
				modal_indices::Array{Int, 1}, point_mat_1D::Array{T, 2};
				scheme = "sparse") where {D, T <: Real}

	cutoff = get_cutoff(scheme, D, n)
	ls::NTuple{D, Int} = ntuple(i->(n+1),D)
	I  = Int[]
	V  = T[]
	c2 = 1

	for l in CartesianRange(ls)
		cutoff(l) && continue

		c2 = inner_loop_to_points(c2, k, l, level, cell, mode, I, V, modal_indices, point_mat_1D)
	end
	return sparsevec(I, V)
end

# Using the above method, builds the transformation matrix
# From the modal sparse basis to
# the "point basis" of values of a given function at the points on the sparse grid
function modal2points(k::Int, n::Int, srefVD::Array{NTuple{3, CartesianIndex{D}}, 1};
	 	scheme = "sparse") where D

	cutoff = get_cutoff(scheme, D, n)
	len = get_size(D, k, n; scheme=scheme)
	modal_ref = D2Vref(1, k, n)
	point_mat_1D = full(modal2points_1D(k, n))
	I = Int[]
	J = Int[]
	V = Real[]
	for c1 in 1:len
		lcm   = srefVD[c1]
		level = lcm[1]
		cell  = lcm[2]
		mode  = lcm[3]
		modal_indices = [modal_ref[(CartesianIndex{1}(level[d],),
								  CartesianIndex{1}(cell[d],),
								  CartesianIndex{1}(mode[d],))] for d in 1:D]
		basis_eval = eval_points(k, n, level, cell, mode,
									modal_indices,
									point_mat_1D;
									scheme=scheme)
		for (c2, val) in zip(findnz(basis_eval)...)
			push!(I, c2)
			push!(J, c1)
			push!(V, val)
		end
	end
	return sparse(I, J, V)
end

# Version of the above function with simpler arguments
function modal2points(D::Int, k::Int, n::Int; scheme="sparse")
	VD = V2Dref(D, k, n; scheme=scheme)
	return modal2points(k, n, VD; scheme=scheme)
end

# Given a specific collocation point, this constructs the corresponding
# combination of nodal basis functions that vanish on all collcation points but that one
function eval_nodal(k::Int, n::Int, level::CartesianIndex{D},
		cell::CartesianIndex{D}, mode::CartesianIndex{D},
		inv_mat_1D::Array{T, 2}; scheme = "sparse") where {D, T <: Real}

	cutoff = get_cutoff(scheme, D, n)
	ls::NTuple{D, Int} = ntuple(i->(n+1),D)
	point_indices = [get_point_index(k, level[d]-1, cell[d], mode[d]) for d in 1:D]

	I = Int[]
	V = T[]
	c2 = 1
	for l in CartesianRange(ls)
		cutoff(l) && continue

		c2 = inner_loop_to_points(c2, k, l, level, cell, mode, I, V, point_indices, inv_mat_1D)
	end
	return sparsevec(I, V)
end

# Using the above, this gives the transformation to go from the point basis
# to the nodal basis of functions
function points2nodal(D::Int, k::Int, n::Int; scheme="sparse")
	cutoff = get_cutoff(scheme, D, n)
	ls	   = ntuple(i->(n+1),D)
	hier_ref = D2Vref(1, k, n)
	inv_mat_1D = full(points2nodal_1D(k, n))
	I = Int[]
	J = Int[]
	V = Real[]

	c1 = 1
	for l in CartesianRange(ls)
		cutoff(l) && continue

		cells = ntuple(i -> 1<<max(0, l[i]-2), D)
		modes = ntuple(i->l[i]==1 ? k : (k-1) ,D)
		for c in CartesianRange(cells)
			for m in CartesianRange(modes)
				coeffs = eval_nodal(k, n, l, c, m,
									inv_mat_1D;
									scheme=scheme)
				for (c2, val) in zip(findnz(coeffs)...)
					push!(I, c2)
					push!(J, c1)
					push!(V, val)
				end
				c1 += 1
			end
		end
	end
	return dropzeros!(sparse(I, J, V))
end


# Given a nodal basis function, this expresses it in terms of the (sparse) DG Basis
function eval_DG(k::Int, n::Int, level::CartesianIndex{D},
		cell::CartesianIndex{D}, mode::CartesianIndex{D},
		hier_ref::Dict{NTuple{3, CartesianIndex{1}}, Int},
		inv_nodal_mat_1D::Array{T, 2}; scheme = "sparse") where {D, T <: Real}

	nodal_indices = [get_point_index(k, level[d]-1, cell[d], mode[d]) for d in 1:D]
	cutoff = get_cutoff(scheme, D, n)
	#len = get_size(D, k, n; scheme=scheme)
	I = Int[]
	V = T[]
	ls::NTuple{D, Int} = ntuple(i->(n+1),D)

	c2 = 1
	for l in CartesianRange(ls)
		cutoff(l) && continue

		c2 = inner_loop_to_modal(c2, k, l, level, cell, mode, I, V, nodal_indices, inv_nodal_mat_1D)
	end
	return sparsevec(I, V)
end

# This gives the transformation matrix between the (sparse) nodal basis
# and the (sparse) DG basis already developed in the previous paper
function nodal2modal(k::Int, n::Int, srefVD::Array{NTuple{3, CartesianIndex{D}}, 1};
		scheme = "sparse") where D

	cutoff = get_cutoff(scheme, D, n)
	ls     = ntuple(i->(n+1), D)
	hier_ref = D2Vref(1, k, n)
	inv_point_mat_1D = full(nodal2modal_1D(k, n))
	I = Int[]
	J = Int[]
	V = Real[]

	c1 = 1
	for l in CartesianRange(ls)
		cutoff(l) && continue

		cells = ntuple(i -> 1<<max(0, l[i]-2), D)
		modes = ntuple(i->l[i]==1 ? k : k-1, D)
		for c in CartesianRange(cells)
			for m in CartesianRange(modes)
				coeffs = eval_DG(k, n, l, c, m,
									hier_ref,
									# srefVD,
									inv_point_mat_1D;
									scheme=scheme)
				for (c2, val) in zip(findnz(coeffs)...)
					push!(I, c2)
					push!(J, c1)
					push!(V, val)
				end
				c1 += 1
			end
		end
	end

	return dropzeros!(sparse(I, J, V))
end

# Version of the above function with simpler arguments
function nodal2modal(D::Int, k::Int, n::Int; scheme="sparse")
	VD = V2Dref(D, k, n; scheme=scheme)
	return nodal2modal(k, n, VD; scheme=scheme)
end

function points2modal(D::Int, k::Int, n::Int)
	return nodal2modal(D, k, n)*points2nodal(D, k, n)
end
