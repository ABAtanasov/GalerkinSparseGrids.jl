# Currently there is an inefficiency in Julia v0.5 and possibly v0.6 with ntuple
# that will hopefully be fixed in the next version. To avoid
# slowdown, we implement a tuple constructor independent of Julia's
# metaprogramming.
#
# The method below is used for building the multi-dimensional
# derivative matrix:
@generated function make_cartesian_index(i::Int, arr1::Int, arr2::CartesianIndex{D}) where D
	levs = [:($q == i ? arr1 : arr2[$q]) for q=1:D]
	quote
		$(Expr(:meta, :inline))
		CartesianIndex{D}($(levs...))
	end
end

# Threshold method for full matrices, producing a sparse one
function threshold(mat::AbstractArray{T, 2},; atol = eps(T)) where T <: Real
	for i in eachindex(mat)
		val = mat[i]
		if abs(val) < atol
			mat[i] = 0
		end
	end
	return sparse(mat)
end

# Threshold method for sparse matrices
function threshold(mat::SparseMatrixCSC{T, Int},; atol = eps(T)) where T <: Real
	(I, J, V) = findnz(mat)
	for (index, val) in enumerate(V)
		if abs(val) < atol
			V[index] = zero(T)
		end
	end
	return dropzeros!(sparse(I, J, V))
end
