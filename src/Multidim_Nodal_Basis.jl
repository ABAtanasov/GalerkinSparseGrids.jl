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


function inner_loop(c2::Int, k::Int, l::CartesianIndex{D}, 
									 level::CartesianIndex{D},
									 cell::CartesianIndex{D},
									 mode::CartesianIndex{D},
									 I::Array{Int, 1},
									 V::Array{T, 1},
									 base_indices::Array{Int, 1},
									 mat_1D::Array{T, 2}) where {D, T<:Real}
	cells::NTuple{D, Int} = ntuple(i -> 1<<max(0, l[i]-2), D)
	modes::NTuple{D, Int} = ntuple(i -> k, D)
	cellrange = CartesianIndices(cells)
	moderange = CartesianIndices(modes)
	for c in cellrange
		!relevant_cell(k, l, c, level, cell) && (c2 += prod(modes); continue)

		for m in moderange
			val = one(T)
			#val = multiply_across(mat_1D, base_indices, k, l, c, m)
			for d in 1:D
				index = get_index_1D(k, l[d]-1, c[d], m[d])
				val *= mat_1D[index, base_indices[d]]
				val == 0 && break
			end
			(abs(val) < eps(T)) && (c2 += 1; continue)
			push!(I, c2)
			push!(V, val)
			c2 += 1
		end
	end
	return c2
end


function make_column(k::Int, n::Int, level::CartesianIndex{D}, 
								cell::CartesianIndex{D}, 
								mode::CartesianIndex{D},
								mat_1D::Array{T, 2};
								scheme="sparse") where {D, T<:Real}
	cutoff = get_cutoff(scheme, D, n)
	ls = ntuple(i->(n+1), D)::NTuple{D, Int}
	base_indices = [get_index_1D(k, level[d]-1, cell[d], mode[d]) for d in 1:D]::Array{Int, 1}

	I = Int[]
	V = T[]
	c2 = 1
	for l in CartesianIndices(ls)
		cutoff(l) && continue

		c2 = inner_loop(c2, k, l, level, cell, mode, I, V, base_indices, mat_1D)
	end
	return sparsevec(I, V)
end


function transform(D::Int, k::Int, n::Int, mat_1D::Array{T,2}; scheme="sparse", atol=1e-14) where {T<:Real}
	cutoff = get_cutoff(scheme, D, n)
	ls	   = ntuple(i->(n+1),D)
	modes  = ntuple(i -> k, D)
	hier_ref = D2Vref(1, k, n)
	I = Int[]
	J = Int[]
	V = Float64[]
	
	c1 = 1
	for l in CartesianIndices(ls)
		cutoff(l) && continue

		cells = ntuple(i -> 1<<max(0, l[i]-2), D)
		for c in CartesianIndices(cells)
			for m in CartesianIndices(modes)
				coeffs = make_column(k, n, l, c, m,  
									mat_1D; 
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
	return threshold(sparse(I, J, V), atol)
end

function transform(D::Int, k::Int, n::Int, from::String, to::String; scheme="sparse", atol=1e-14)
	if Set{String}([from, to]) == Set{String}(["modal", "points"])
		T1 = transform(D, k, n, from, "nodal"; scheme=scheme, atol=atol)
		T2 = transform(D, k, n, "nodal", to; scheme=scheme, atol=atol)
		return threshold(T2 * T1, atol)
	else
		mat_1D = Matrix(transform_1D(k, n, from, to))
		return transform(D, k, n, mat_1D; scheme=scheme, atol=atol)
	end
end

