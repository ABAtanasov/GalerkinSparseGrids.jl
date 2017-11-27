using GalerkinSparseGrids
using Base.Test
using Cubature

#---------------------------------------
# Testing vector hier DG reconstruction
#---------------------------------------

print("Testing full DG reconstruction 1-D with vector coefficients... ")

for k in 1:5
    for l in 1:6
        vect = vcoeffs_DG(1, k, l, x->sin(4*x[1]); scheme="full")
        dict = V2D(1, k, l, vect; scheme="full")
        err = x->(reconstruct_DG(dict,[x[1]])-sin(4*x[1]))^2
        @test hquadrature(err,0,1; reltol=1.0e-9, abstol=1.0e-12)[1] < 1/(1<<((l+k-1)))
    end
end

println("Test Passed.")

# Testing D2V and V2D

# Full case:

print("Testing D2V and V2D Full Case 2D... ")

for k in 1:5
    for l in 1:3
        dict = coeffs_DG(2, k, l, x->sin(4*x[1]+x[2]); scheme="full")
        vect = D2V(2, k, l, dict; scheme="full")
        @test V2D(2, k, l, vect; scheme="full")==dict
    end
end

for k in 1:5
    for l in 1:3
        vect = vcoeffs_DG(2, k, l, x->sin(4*x[1]+x[2]); scheme="full")
        dict = V2D(2, k, l, vect; scheme="full")
        @test D2V(2, k, l, dict; scheme="full")==vect
    end
end

println("Test Passed.")

# Sparse case:

print("Testing D2V and V2D Sparse Case 2D... ")

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

println("Test Passed.")
