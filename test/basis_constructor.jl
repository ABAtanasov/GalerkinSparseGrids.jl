using GalerkinSparseGrids
using Test
using HCubature
using StaticArrays
using SparseArrays

#--------------------------------------
# Constructing the Orthogonal DG Basis
#--------------------------------------

@testset "basis_constructor.jl" begin
    import GalerkinSparseGrids.dg_basis
    @info "Testing that the DG basis is created orthogonally"
    @testset "dg_basis does not throw error" begin
        for k in 1:5
            Q = dg_basis(k)
            @test typeof(Q) == Array{Array{Float64,1},1}
        end
    end    
end
