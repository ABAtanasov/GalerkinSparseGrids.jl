using GalerkinSparseGrids
using Base.Test
using Cubature

#---------------------------------------
# Testing regular hier DG reconstruction
#---------------------------------------

print("Testing full discontinuous Galerkin (DG) reconstruction 1-D... ")

for k in 1:5
    for l in 1:6
        dict = coeffs_DG(1, k, l, x->sin(4*x[1]); scheme="full")
        err = x->(reconstruct_DG(dict, [x])-sin(4*x))^2
        @test hquadrature(err, 0, 1; abstol=1.0e-10)[1] < 1/(1<<(l+k-1))
    end
end

for k in 1:5
    for l in 1:6
        dict = coeffs_DG(1, k, l, x->sin(4*x[1]); scheme="sparse")
        err = x->(reconstruct_DG(dict, [x])-sin(4*x))^2
        @test hquadrature(err, 0, 1; abstol=1.0e-10)[1] < 1/(1<<(l+k-1))
    end
end

println("Test Passed.")

# Sparsification should be trivial in 1D

print("Testing that both reconstructions are equivalent... ")

for k in 1:5
    for l in 1:6
        full_coeffs = coeffs_DG(1, k, l, x->sin(4*x[1]); scheme="full")
        sparse_coeffs = coeffs_DG(1, k, l, x->sin(4*x[1]); scheme="full")
        diff_full = x->(reconstruct_DG(full_coeffs, [x])-sin(4*x))^2
        diff_sparse = x->(reconstruct_DG(sparse_coeffs, [x])-sin(4*x))^2
        diff_both = x->(reconstruct_DG(full_coeffs, [x])-reconstruct_DG(sparse_coeffs, [x]))^2

        @test abs(hquadrature(diff_full,0,1;abstol=1.0e-9)[1]-
                  hquadrature(diff_sparse,0,1;abstol=1.0e-9)[1]) < 1.0e-14
        @test hquadrature(diff_both,0,1;abstol=1.0e-9)[1] < 1.0e-14
    end
end

println("Test Passed.")

# Testing sparse DG reconstruction in multidimensional space:

print("Testing sparse DG basis reconstruction 2-D... ")

for k in 1:5
    for l in 1:6
        dict = coeffs_DG(2, k, l, x->sin(4*x[1]+x[2]))
        err = x->(reconstruct_DG(dict, [x[1],x[2]])-sin(4*x[1]+x[2]))^2
        @test hcubature(err, [0,0], [1,1]; abstol=1.0e-10, maxevals=500)[1] < 1/(1<<(l+k-2))
    end
end

println("Test Passed.")
