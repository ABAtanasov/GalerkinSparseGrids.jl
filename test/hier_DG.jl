using GalerkinSparseGrids
using Test
using HCubature
using StaticArrays
using SparseArrays


#---------------------------------------
# Testing regular hier DG reconstruction
#---------------------------------------

@info "Testing full discontinuous Galerkin (DG) reconstruction 1-D... "
@testset "hier_DG.jl" begin
    @testset "Regular hier DG reconstruction" begin
        for k in 1:5
            for l in 1:6
                dict = coeffs_DG(1, k, l, x->sin(4*x[1]); scheme="full")
                err = x->(reconstruct_DG(dict, [x])-sin(4*x))^2
                @test hquadrature(err, 0, 1; atol=1.0e-10)[1] < 1/(1<<(l+k-1))
            end
        end

        for k in 1:5
            for l in 1:6
                dict = coeffs_DG(1, k, l, x->sin(4*x[1]); scheme="sparse")
                err = x->(reconstruct_DG(dict, [x])-sin(4*x))^2
                @test hquadrature(err, 0, 1; atol=1.0e-10)[1] < 1/(1<<(l+k-1))
            end
        end
    end

    # Sparsification should be trivial in 1D
    @info "Testing that both reconstructions are equivalent... "
    @testset "Equivalent Reconstructions" begin
        for k in 1:5
            for l in 1:6
                full_coeffs = coeffs_DG(1, k, l, x->sin(4*x[1]); scheme="full")
                sparse_coeffs = coeffs_DG(1, k, l, x->sin(4*x[1]); scheme="full")
                diff_full = x->(reconstruct_DG(full_coeffs, [x])-sin(4*x))^2
                diff_sparse = x->(reconstruct_DG(sparse_coeffs, [x])-sin(4*x))^2
                diff_both = x->(reconstruct_DG(full_coeffs, [x])-reconstruct_DG(sparse_coeffs, [x]))^2

                @test abs(hquadrature(diff_full,0,1; atol=1.0e-9)[1]-
                          hquadrature(diff_sparse,0,1; atol=1.0e-9)[1]) < 1.0e-14
                @test hquadrature(diff_both,0,1; atol=1.0e-9)[1] < 1.0e-14
            end
        end
    end

# Testing sparse DG reconstruction in multidimensional space:

    @info "Testing sparse DG basis reconstruction 2-D..."
    @testset "Sparse DG reconstruction in multidimensional space" begin
        for k in 1:5
            for l in 1:6
                dict = coeffs_DG(2, k, l, x->sin(4*x[1]+x[2]))
                err = x->(reconstruct_DG(dict, [x[1],x[2]])-sin(4*x[1]+x[2]))^2
                @test hcubature(err, [0,0], [1,1]; atol=1.0e-10, maxevals=500)[1] < 1/(1<<(l+k-2))
            end
        end
    end
end
