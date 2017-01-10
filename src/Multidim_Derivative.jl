#------------------------------------------------------------
#
# Construction of multidimensional DG derivative matrices
#
#------------------------------------------------------------

# Efficiency criticality: LOW
# Computations only performed once

# Accuracy criticality: HIGH
# Critical for accurate PDE evolution

function full_D_matrix{D}(i::Int, k::Int, n::Int,
    srefVD::Array{CartesianIndex{D},2}, srefDV::Dict{Array{CartesianIndex{D},1},Int})
    
    V2D_1D = full_referenceV2D(k,(n+1,))
    D2V_1D = full_referenceD2V(k,(n+1,))
    Dmat_1D = periodic_hier_DLF_Matrix(0, k, n)

    len = length(srefVD[:,1])
    sMat= spzeros(len,len)

    for c1 in 1:len
        lpf = view(srefVD,c1,:)
        l = lpf[1][i]
        p = lpf[2][i]
        f = lpf[3][i]
        vc1 = D2V_1D[[CartesianIndex((l,)),CartesianIndex((p,)),CartesianIndex((f,))]]
        dvc1s = view(Dmat_1D,:, vc1)
        
        for j in 1:length(dvc1s)
            if abs(dvc1s[j])<=1.0e-15
                continue
            end

            lpf2 = view(V2D_1D,j,:)
            
            level2 = CartesianIndex{D}(ntuple(j-> j==i?lpf2[1][1]:lpf[1][j] ,D))
            place2 = CartesianIndex{D}(ntuple(j-> j==i?lpf2[2][1]:lpf[2][j] ,D))
            f_number2 = CartesianIndex{D}(ntuple(j-> j==i?lpf2[3][1]:lpf[3][j] ,D))
            
            c2 = srefDV[[level2, place2, f_number2]]
            sMat[c2,c1] += dvc1s[j]
        end
    end
    return sMat
end

function sparse_D_matrix{D}(i::Int, k::Int, n::Int,
    srefVD::Array{CartesianIndex{D},2}, srefDV::Dict{Array{CartesianIndex{D},1},Int})
    V2D_1D = full_referenceV2D(k,(n+1,))
    D2V_1D = full_referenceD2V(k,(n+1,))
    Dmat_1D = periodic_hier_DLF_Matrix(0, k, n)
    
    len = length(srefVD[:,1])
    sMat= spzeros(len,len)
    for c1 in 1:len
        lpf = view(srefVD,c1,:)
        l = lpf[1][i]
        p = lpf[2][i]
        f = lpf[3][i]
        vc1 = D2V_1D[[CartesianIndex((l,)),CartesianIndex((p,)),CartesianIndex((f,))]]
        dvc1s = view(Dmat_1D,:, vc1)
        
        for j in 1:length(dvc1s)
            if abs(dvc1s[j])<=1.0e-15
                continue
            end
            lpf2 = view(V2D_1D,j,:)
            level2 = CartesianIndex{D}(ntuple(j-> j==i?lpf2[1][1]:lpf[1][j] ,D))
            diag_level=0;
            for q in 1:D
                diag_level+=level2[q]
            end
            if diag_level >= n + D 
                continue
            end
            place2 = CartesianIndex{D}(ntuple(j-> j==i?lpf2[2][1]:lpf[2][j] ,D))
            f_number2 = CartesianIndex{D}(ntuple(j-> j==i?lpf2[3][1]:lpf[3][j] ,D))
            
            c2 = srefDV[[level2, place2, f_number2]]
            sMat[c2,c1] += dvc1s[j]
        end
    end
    return sMat
end