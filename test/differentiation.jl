using GalerkinSparseGrids
using Test
using HCubature
using SparseArrays

#--------------------------------------
# Testing Differentiation
#--------------------------------------

@testset "differentiation.jl" begin 
    @testset "Testing differentiation 1-D DG basis... " begin
        k=3
        for l in 2:5
            frefVD = V2Dref(1, k, l);
            frefDV = D2Vref(1, k, l);
            D_op = D_matrix(1, k, l, frefVD, frefDV; scheme="full")
            vcoeffs = vcoeffs_DG(1, k, l, x->cos(2*pi*x[1]); scheme="full")
            dvcoeffs = *(D_op,vcoeffs)
            dict= V2D(1, k, l, vcoeffs; scheme="full")
            ddict = V2D(1, k, l, dvcoeffs; scheme="full")
            err = x->(reconstruct_DG(ddict,[x[1]])+2*pi*sin(2*pi*x[1]))^2
            @test hquadrature(err, 0, 1, abstol=1.0e-10, maxevals=500)[1] < 1/(1<<(k+l-2))
        end
    end

    @testset "Testing differentiation 2-D full DG basis... " begin
        D = 2
        k = 3
        for l in 2:5
            D_op = D_matrix(D, 1, k, l; scheme="full")
            vcoeffs = vcoeffs_DG(D, k, l, x->cos(2*pi*x[1])*cos(2*pi*x[2]); scheme="full")
            dvcoeffs = *(D_op, vcoeffs)
            dict= V2D(D, k, l, vcoeffs; scheme="full")
            ddict = V2D(D, k, l, dvcoeffs; scheme="full")
            err = x->(reconstruct_DG(ddict,[x[1],x[2]])+2*pi*sin(2*pi*x[1])*cos(2*pi*x[2]))^2
            @test hcubature(err, [0,0], [1,1], abstol=1.0e-10, maxevals=500)[1] < 1/(1<<(k+l-2))
        end
    end
end
