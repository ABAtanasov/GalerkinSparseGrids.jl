#------------------------------------------------------------
#
# Construction of multidimensional DG derivative matrices
#
#------------------------------------------------------------

# Efficiency criticality: HIGH
# Matrix computations are only performed once, 
# but this is currently the bottleneck

# Accuracy criticality: HIGH
# Critical for accurate PDE evolution


@generated function make_cartesian_index{D}(i::Int, arr1::Int, arr2::CartesianIndex{D})
	levs = [:($q == i ? arr1 : arr2[$q]) for q=1:D]
	quote
		$(Expr(:meta, :inline))
		CartesianIndex{D}($(levs...))
	end
end

function sparse_D_matrix{D}(i::Int, k::Int, n::Int,
    						srefVD::Array{NTuple{3,CartesianIndex{D}},1}, 
							srefDV::Dict{NTuple{3,CartesianIndex{D}},Int})
	
    V2D_1D = full_referenceV2D(k,(n+1,))
    D2V_1D = full_referenceD2V(k,(n+1,))
    Dmat_1D = periodic_hier_DLF_Matrix(0, k, n)
    
    len = length(srefVD[:,1])
	I = Int[]
	J = Int[]
	K = Float64[]
	
    for c1 in 1:len
        lpf = srefVD[c1] #view(srefVD,c1,:)
        l = lpf[1][i]
        p = lpf[2][i]
        f = lpf[3][i]
        vc1 = D2V_1D[(CartesianIndex(l),CartesianIndex(p),CartesianIndex(f))]
        dvc1s = view(Dmat_1D,:, vc1)
        
        for j in 1:length(dvc1s)
			lpf2 = V2D_1D[j]  #view(V2D_1D,j,:)
			
			# Currently there is an inefficiency in Julia v0.5 with ntuple
			# that will hopefully be fixed in the next version. To avoid
			# slowdown, we implement a tuple constructor independent of Julia's
			# metaprogramming:
			
            level2 = make_cartesian_index(i, lpf2[1][1], lpf[1])
            diag_level=0
            for q in 1:D
                diag_level += level2[q]
            end
            if diag_level >= n + D 
                continue
            end
            place2 = make_cartesian_index(i, lpf2[2][1], lpf[2])
			f_number2 = make_cartesian_index(i, lpf2[3][1], lpf[3])
			#place2 = CartesianIndex{D}(ntuple(j-> j==i?lpf2[2][1]:lpf[2][j] ,D))
            #f_number2 = #CartesianIndex{D}(ntuple(j-> j==i?lpf2[3][1]:lpf[3][j] ,D))
            
            c2 = srefDV[(level2, place2, f_number2)]
			push!(I, c2)
			push!(J, c1)
			push!(K, dvc1s[j])
        end
    end
	
	Dfinal = sparse(I, J ,K, len, len, +)
	# Does not seem helpful for this matrix:
	# dropzeros!(Dfinal)
	
    return Dfinal
end


function full_D_matrix{D}(i::Int, k::Int, n::Int,
    						frefVD::Array{NTuple{3,CartesianIndex{D}},1}, 
							frefDV::Dict{NTuple{3,CartesianIndex{D}},Int})
	
    V2D_1D = full_referenceV2D(k,(n+1,))
    D2V_1D = full_referenceD2V(k,(n+1,))
    Dmat_1D = periodic_hier_DLF_Matrix(0, k, n)
    
    len = length(frefVD[:,1])
	I = Int[]
	J = Int[]
	K = Float64[]
	
    for c1 in 1:len
        lpf = frefVD[c1]
        l = lpf[1][i]
        p = lpf[2][i]
        f = lpf[3][i]
        vc1 = D2V_1D[(CartesianIndex(l),CartesianIndex(p),CartesianIndex(f))]
        dvc1s = view(Dmat_1D,:, vc1)
        
        for j in 1:length(dvc1s)
			lpf2 = V2D_1D[j]  #view(V2D_1D,j,:)

            level2 = make_cartesian_index(i, lpf2[1][1], lpf[1])
            place2 = make_cartesian_index(i, lpf2[2][1], lpf[2])
			f_number2 = make_cartesian_index(i, lpf2[3][1], lpf[3])
            
            c2 = frefDV[(level2, place2, f_number2)]
			push!(I, c2)
			push!(J, c1)
			push!(K, dvc1s[j])
        end
    end
	
	Dfinal = sparse(I, J ,K, len, len, +)
		
    return Dfinal
end
