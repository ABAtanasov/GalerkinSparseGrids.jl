#------------------------------------------------------------
#
# Construction of multidimensional DG derivative matrices
# based on the (precomputed) 1-D derivative matrix
#
#------------------------------------------------------------

# Efficiency criticality: HIGH
# Matrix computations are only performed once, 
# but this is currently the bottleneck

# Accuracy criticality: HIGH
# Critical for accurate PDE evolution

function D_matrix{D}(i::Int, k::Int, n::Int,
							srefVD::Array{NTuple{3,CartesianIndex{D}},1},
							srefDV::Dict{NTuple{3,CartesianIndex{D}},Int};
							scheme="sparse")
	
	cutoff = get_cutoff(scheme, D, n)
	
	V2D_1D = V2Dref(1,k,n)
	D2V_1D = D2Vref(1,k,n)
	Dmat_1D = periodic_DLF_Matrix(k, n)
	
	len = length(srefVD[:,1])
	I = Int[]
	J = Int[]
	V = Float64[]
	
	for c1 in 1:len
		lpf = srefVD[c1]
		l = lpf[1][i]
		p = lpf[2][i]
		f = lpf[3][i]
		vc1 = D2V_1D[(CartesianIndex(l),CartesianIndex(p),CartesianIndex(f))]
		dvc1s = view(Dmat_1D,:, vc1)
		
		for j in 1:length(dvc1s)
			lpf2 = V2D_1D[j]
			# using tuple constructor independent of
			# Julia's metaprogramming:
			level2 = make_cartesian_index(i, lpf2[1][1], lpf[1])
			cutoff(level2) && continue
			
			place2 = make_cartesian_index(i, lpf2[2][1], lpf[2])
			f_number2 = make_cartesian_index(i, lpf2[3][1], lpf[3])
			c2 = srefDV[(level2, place2, f_number2)]
			push!(I, c2)
			push!(J, c1)
			push!(V, dvc1s[j])
		end
	end
	# dropzeros! does not seem helpful for this matrix:
	return sparse(I, J, V, len, len, +)
end


function D_matrix(D::Int, i::Int, k::Int, n::Int; scheme="sparse")
	VD = V2Dref(D, k, n; scheme=scheme)
	DV = D2Vref(D, k, n; scheme=scheme)
	return D_matrix(i, k, n, VD, DV; scheme=scheme)
end

function grad_matrix(D::Int, k::Int, n::Int; scheme="sparse")
	return [D_matrix(D, i, k, n; scheme=scheme) for i in 1:D]
end

function laplacian_matrix(D::Int, k::Int, n::Int; scheme="sparse")
	len = get_size(D, k, n; scheme=scheme)	
	laplac = spzeros(len, len)
	for i in 1:D
		D_op = D_matrix(D, i, k, n; scheme=scheme)
		laplac += *(D_op,D_op)
	end
	return laplac
end