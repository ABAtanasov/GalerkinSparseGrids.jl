using GalerkinSparseGrids
using Base.Test
using Cubature

#--------------------------------------
# Testing hat reconstruction
#--------------------------------------

# Position:
print("Testing position hat basis reconstruction 1-D... ")

for l in 1:7
    pos_coeffs = standard_coeffs(x->sin(4*x[1]),(l,))
    @test hquadrature(x->(standard_reconstruct(pos_coeffs, (l,), (x,))-sin(4*x))^2,0,1;abstol=1.0e-9)[1]<1/(1<<(l))
end

println("Test Passed.")

# Hierarchical:

print("Testing full hat basis reconstruction 1-D... ")

for l in 3:9
    hier_coeffs = coeffs_hat(1, l, x->sin(4*x[1]))
    @test hquadrature(x->(reconstruct_hat(hier_coeffs, (x,))-sin(4*x))^2,0,1;abstol=1.0e-9)[1]<1/(1<<(l))
end

println("Test Passed.")

# Hierarchical and position basis should yield the same results for hat functions:

print("Testing that both reconstructions are equivalent... ")

for l in 1:7
    pos_coeffs=standard_coeffs(x->sin(4*x[1]),(l,))
    hier_coeffs = coeffs_hat(1, l+1, x->sin(4*x[1]))
    diff_standard = x->(standard_reconstruct(pos_coeffs, (l,), (x,))-sin(4*x))^2
    diff_hier = x->(reconstruct_hat(hier_coeffs, (x,))-sin(4*x))^2
    diff_both = x->(standard_reconstruct(pos_coeffs, (l,), (x,))-reconstruct_hat(hier_coeffs, (x,)))^2

    @test abs(hquadrature(diff_standard,0,1;abstol=1.0e-9)[1]-
              hquadrature(diff_hier,0,1;abstol=1.0e-9)[1])<1.0e-14
    @test hquadrature(diff_both,0,1;abstol=1.0e-9)[1]<1.0e-14
end

println("Test Passed.")

# Multidimensional hat reconstruction:

print("Testing position hat basis reconstruction 2-D... ")

for l in 1:5

    pos_coeffs = standard_coeffs(x->sin(4*x[1]+x[2]),(l,l))
    pos_err = x->(standard_reconstruct(pos_coeffs, (l,l), (x[1],x[2]))-sin(4*x[1]+x[2]))^2
    @test hcubature(pos_err, [0,0], [1,1]; abstol=1.0e-9, maxevals=500)[1]<1/(1<<(l))
end

println("Test Passed.")

print("Testing full hat basis reconstruction 2-D... ")

for l in 3:7
    hier_coeffs = coeffs_hat(2, l, x->sin(4*x[1]+x[2]); scheme="full")
    hier_err = x->(reconstruct_hat(hier_coeffs, (x[1],x[2]))-sin(4*x[1]+x[2]))^2
    @test hcubature(hier_err, [0,0], [1,1]; abstol=1.0e-9)[1]<1/(1<<(l))
end

println("Test Passed.")

print("Testing sparse hat basis reconstruction 2-D... ")

for l in 1:5
    sparse_coeffs = coeffs_hat(2, l, x->sin(4*x[1]+x[2]); scheme="sparse")
    sparse_err = x->(reconstruct_hat(sparse_coeffs, (x[1],x[2]))-sin(4*x[1]+x[2]))^2
    @test hcubature(sparse_err, [0,0], [1,1];abstol=1.0e-9)[1]<1/(1<<(l))
end

println("Test Passed.")
