using GalerkinSparseGrids
using Test
using HCubature
using StaticArrays
using SparseArrays

#---------------------------------------
# Testing vector hier DG reconstruction
#---------------------------------------

@info "Testing full DG reconstruction 1-D with vector coefficients... "

@testset "vhier_DC.jl" begin
	@info "Testing full DG reconstruction with vector coefficients in 1D"
    @testset "Full DG reconstruction 1D with vector coefficients" begin
        for k in 1:5
            for l in 1:6
                vect = vcoeffs_DG(1, k, l, x->sin(4*x[1]); scheme="full")
                dict = V2D(1, k, l, vect; scheme="full")
                err = x->(reconstruct_DG(dict,[x[1]])-sin(4*x[1]))^2
                @test hquadrature(err,0,1; rtol=1.0e-9, atol=1.0e-12)[1] < 1/(1<<((l+k-1)))
            end
        end
    end

    # Testing D2V and V2D

    # Full case:

    @info "Testing D2V and V2D Full Case 2D... "
    @testset "D2V and V2D full case 2D" begin
        for k in 1:5
            for l in 1:3
                dict = coeffs_DG(2, k, l, x->sin(4*x[1]+x[2]); scheme="full")
                vect = D2V(2, k, l, dict; scheme="full")
                @test V2D(2, k, l, vect; scheme="full") == dict
            end
        end

        for k in 1:5
            for l in 1:3
                vect = vcoeffs_DG(2, k, l, x->sin(4*x[1]+x[2]); scheme="full")
                dict = V2D(2, k, l, vect; scheme="full")
                @test D2V(2, k, l, dict; scheme="full") == vect
            end
        end
    end

    # Sparse case:

    @info "Testing D2V and V2D Sparse Case 2D... "
    @testset "D2V and V2D sparse case 2D" begin
        for k in 1:5
            for l in 1:3
                dict = coeffs_DG(2, k, l, x->sin(4*x[1]+x[2]); scheme="sparse")
                vect = D2V(2, k, l, dict; scheme="sparse")
                @test V2D(2, k, l, vect; scheme="sparse")==dict
            end
        end

        for k in 1:5
            for l in 1:3
                vect = vcoeffs_DG(2, k, l, x->sin(4*x[1]+x[2]); scheme="sparse")
                dict= V2D(2, k, l, vect; scheme="sparse")
                @test D2V(2, k, l, dict; scheme="sparse")==vect
            end
        end
    end
end
