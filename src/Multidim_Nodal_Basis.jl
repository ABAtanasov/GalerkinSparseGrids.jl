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
							   cell1 ::CartesianIndex{D}, 
							   level2::CartesianIndex{D}, 
							   cell2 ::CartesianIndex{D}) where D
	bool = true
	for d in 1:D
		bool &= relevant_cell_1D(k, level1[d]-1, cell1[d], level2[d]-1, cell2[d])
	end
	return bool
end


function inner_loop(i::Int, j::Int, I::Array{Int, 1}, J::Array{Int, 1}, V::Array{Float64, 1},
					k::Int, l2::CartesianIndex{D}, 
					l1::CartesianIndex{D}, c1::CartesianIndex{D}, m1::CartesianIndex{D},
					base_indices::Array{Int, 1}, mat_1D::Array{T, 2}) where {D, T<:Real}
	cells::NTuple{D, Int} = ntuple(i -> 1<<max(0, l2[i]-2), D)
	modes::NTuple{D, Int} = ntuple(i -> k, D)
	cellrange = CartesianIndices(cells)
	moderange = CartesianIndices(modes)
	for c2 in cellrange
		!relevant_cell(k, l2, c2, l1, c1) && (i += prod(modes); continue)

		for m2 in moderange
			val = one(T)
			#val = multiply_across(mat_1D, base_indices, k, l, c, m)
			for d in 1:D
				index = get_index_1D(k, l2[d]-1, c2[d], m2[d])
				val *= mat_1D[index, base_indices[d]]
				val == 0 && break
			end
			(abs(val) < eps(T)) && (i += 1; continue)
			push!(I, i)
			push!(J, j)
			push!(V, val)
			i += 1
		end
	end
	return i
end


function make_column(j::Int, I::Array{Int, 1}, J::Array{Int, 1}, V::Array{Float64, 1},
					 k::Int, n::Int, 
					 l1::CartesianIndex{D}, c1::CartesianIndex{D}, m1::CartesianIndex{D},
					 mat_1D::Array{T, 2}; scheme="sparse") where {D, T<:Real}
	cutoff = get_cutoff(scheme, D, n)
	ls = ntuple(i -> (n+1), D)::NTuple{D, Int}
	base_indices = [get_index_1D(k, l1[d]-1, c1[d], m1[d]) for d in 1:D]::Array{Int, 1}

	i = 1
	for l2 in CartesianIndices(ls)
		cutoff(l2) && continue
		i = inner_loop(i, j, I, J, V, k, l2, l1, c1, m1, base_indices, mat_1D)
	end
end


function transform(D::Int, k::Int, n::Int, mat_1D::Array{T,2}; scheme="sparse", atol=1e-12) where {T<:Real}
	cutoff = get_cutoff(scheme, D, n)
	levels = ntuple(i -> (n+1),D)
	modes  = ntuple(i -> k, D)
	I = Int[]
	J = Int[]
	V = Float64[]
	
	j = 1
	for l1 in CartesianIndices(levels)
		cutoff(l1) && continue

		cells1 = ntuple(i -> 1<<max(0, l1[i]-2), D)
		for c1 in CartesianIndices(cells1)
			for m1 in CartesianIndices(modes)
				coeffs = make_column(j, I, J, V, k, n, l1, c1, m1, mat_1D; scheme=scheme)
				j += 1
			end
		end
	end
	return threshold(sparse(I, J, V), atol)
end

function transform(D::Int, k::Int, n::Int, from::String, to::String; scheme="sparse", atol=1e-12)
	# Its better not to multiply the matrices from modal -> nodal with nodal -> points 
	# or vice-versa. Instead, have both transform matrices in memory
	# and multiply a vector once by each of the two to get to the other basis - 
	# So this if-statement is not an efficient thing to do in general
	if Set{String}([from, to]) == Set{String}(["modal", "points"])
		T1 = transform(D, k, n, from, "nodal"; scheme=scheme, atol=atol)
		T2 = transform(D, k, n, "nodal", to; scheme=scheme, atol=atol)
		return threshold(T2 * T1, atol)
	else
		mat_1D = Matrix(transform_1D(k, n, from, to))
		return transform(D, k, n, mat_1D; scheme=scheme, atol=atol)
	end
end


# --------------------------------------------------
# Below this line is for comparison testing only -
# Not for use 
# --------------------------------------------------

function transform2(D::Int, k::Int, n::Int, mat_1D::Array{T,2}; scheme="sparse", atol=1e-12) where {T<:Real}
	cutoff = get_cutoff(scheme, D, n)
	levels = ntuple(i->(n+1),D)
	modes  = ntuple(i -> k, D)
	hier_ref = D2Vref(1, k, n)
	I = Int[]; J = Int[]; V = Float64[]
	i::Int = 1; j::Int = 1; val::Float64 = one(Float64)
	index::Int = 1
	
	levelrange = CartesianIndices(levels); moderange = CartesianIndices(modes)
	base_indices::Array{Int,1} = zeros(Int, D)
	cells1::Array{Int,1} = zeros(Int, D); 
	cells2::Array{Int,1} = zeros(Int, D)
	for l1 in levelrange
		cutoff(l1) && continue

		cells1 = [1<<max(0, l1[q]-2) for q in 1:D]
		for c1 in CartesianIndices((cells1...,))
			for m1 in moderange
				#1D indices corresponding to each (l1[d], c1[d], m1[d])
				base_indices = [get_index_1D(k, l1[d]-1, c1[d], m1[d]) for d in 1:D]

				i = 1
				for l2 in levelrange
					cutoff(l2) && continue

					cells2 = [1<<max(0, l2[q]-2) for q in 1:D]
					for c2 in CartesianIndices((cells2...,))
						# If there's no overlap - skip
						!relevant_cell(k, l2, c2, l1, c1) && (i += prod(modes); continue)
						for m2 in moderange
							val = one(Float64)
							for d in 1:D
								index = get_index_1D(k, l2[d]-1, c2[d], m2[d])
								val *= mat_1D[index, base_indices[d]]
								val == 0 && break
							end
							(abs(val) < eps(T)) && (i += 1; continue)
							push!(I, i)
							push!(J, j)
							push!(V, val)
							i += 1
						end
					end
				end
				j += 1
			end
		end
	end
	return threshold(sparse(I, J, V), atol)
end

function transform2(D::Int, k::Int, n::Int, from::String, to::String; scheme="sparse", atol=1e-12)
	if Set{String}([from, to]) == Set{String}(["modal", "points"])
		T1 = transform2(D, k, n, from, "nodal"; scheme=scheme, atol=atol)
		T2 = transform2(D, k, n, "nodal", to; scheme=scheme, atol=atol)
		return threshold(T2 * T1, atol)
	else
		mat_1D = Matrix(transform_1D(k, n, from, to))
		return transform2(D, k, n, mat_1D; scheme=scheme, atol=atol)
	end
end
