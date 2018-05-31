# Currently there is an inefficiency in Julia v0.5 and possibly v0.6 with ntuple
# that will hopefully be fixed in the next version. To avoid
# slowdown, we implement a tuple constructor independent of Julia's
# metaprogramming. 
#
# The method below is used for building the multi-dimensional 
# derivative matrix:
@generated function make_cartesian_index{D}(i::Int, arr1::Int, arr2::CartesianIndex{D})
	levs = [:($q == i ? arr1 : arr2[$q]) for q=1:D]
	quote
		$(Expr(:meta, :inline))
		CartesianIndex{D}($(levs...))
	end
end

# Threshold method for full matrices, producing a sparse one
function threshold{T<:Real}(mat::Array{T, 2}; abs_tol=eps(T))
	for i in eachindex(mat)
		val = mat[i]
		if abs(val) < abs_tol
			mat[i] = 0
		end
	end
	return sparse(mat)
end

# Threshold method for sparse matrices
function threshold{T<:Real}(mat::SparseMatrixCSC{T, Int}; abs_tol=eps(T))
	(I, J, V) = findnz(mat)
	for (index, val) in enumerate(V)
		if abs(val) < abs_tol
			V[index] = zero(T)
		end
	end
	return dropzeros!(sparse(I, J, V))
end