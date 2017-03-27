using GalerkinSparseGrids
using Base.Test
using Cubature
using ODE

println("Beginning Tests of GalerkinSparseGrids.jl--------------------------------------")

# #--------------------------------------
# # Elementary Tests
# #--------------------------------------
#
# import GalerkinSparseGrids.pos
# print("Testing pos function... ")
# for i in -10:0
# 	@test pos(i)==0
# end
# for i in 0:10
# 	@test pos(i)==i
# end
# println("Tests Passed.")
#
# print("Testing hat_index function... ")
# import GalerkinSparseGrids.hat_index
# for l in 1:5
# 	@test hat_index(1.3,l)==(1<<(l-1))
# end
# for l in 1:5
# 	@test hat_index(.01,l)==1
# end
# @test hat_index(.3,1)==1
# @test hat_index(.3,2)==1
# @test hat_index(.3,3)==2
# @test hat_index(.3,4)==3
# @test hat_index(.3,5)==5
# println("Tests Passed.")
#
# import GalerkinSparseGrids.inner_product
# print("Testing 1D inner product... ")
# @test abs(inner_product(x->x[1]^2, x->x[1]^3, (0,),CartesianIndex((1,)))-(1/6)) < 2*eps(Float64)
# @test abs(inner_product(x->sin(pi*x[1]), x->cos(pi*x[1]), (0,),CartesianIndex((1,)))) < eps(Float64)
# println("Tests Passed.")
# print("Testing 2D inner product... ")
# @test abs(inner_product(x->(x[1]^2+x[2]^2), x->x[1]^3, (0,0),CartesianIndex((1,1)))-.25)<2*eps(Float64)
# println("Tests Passed.")
# print("Testing 3D inner product... ")
# @test abs(inner_product(x->(x[1]^2+x[2]^2-x[3]^2), x->x[1]^3+x[2]+x[3], (0,0,0),CartesianIndex((1,1,1)))-.5)<2*eps(Float64)
# println("Tests Passed.")
#
# #--------------------------------------
# # Testing hat reconstruction
# #--------------------------------------
#
# # Position:
# print("Testing position hat basis reconstruction 1-D... ")
#
# for l in 1:7
#     pos_coeffs = standard_coeffs(x->sin(4*x[1]),(l,))
#     @test hquadrature(x->(standard_reconstruct(pos_coeffs, (l,), (x,))-sin(4*x))^2,0,1;abstol=1.0e-9)[1]<1/(1<<(l))
# end
#
# println("Test Passed.")
#
# # Hierarchical:
#
# print("Testing full hat basis reconstruction 1-D... ")
#
# for l in 3:9
#     hier_coeffs = coeffs_hat(1, l, x->sin(4*x[1]))
#     @test hquadrature(x->(reconstruct_hat(hier_coeffs, (x,))-sin(4*x))^2,0,1;abstol=1.0e-9)[1]<1/(1<<(l))
# end
#
# println("Test Passed.")
#
# # Hierarchical and position basis should yield the same results for hat functions:
#
# print("Testing that both reconstructions are equivalent... ")
#
# for l in 1:7
#     pos_coeffs=standard_coeffs(x->sin(4*x[1]),(l,))
#     hier_coeffs = coeffs_hat(1, l+1, x->sin(4*x[1]))
#     diff_standard = x->(standard_reconstruct(pos_coeffs, (l,), (x,))-sin(4*x))^2
#     diff_hier = x->(reconstruct_hat(hier_coeffs, (x,))-sin(4*x))^2
#     diff_both = x->(standard_reconstruct(pos_coeffs, (l,), (x,))-reconstruct_hat(hier_coeffs, (x,)))^2
#
#     @test abs(hquadrature(diff_standard,0,1;abstol=1.0e-9)[1]-
#               hquadrature(diff_hier,0,1;abstol=1.0e-9)[1])<1.0e-14
#     @test hquadrature(diff_both,0,1;abstol=1.0e-9)[1]<1.0e-14
# end
#
# println("Test Passed.")
#
# # Multidimensional hat reconstruction:
#
# print("Testing position hat basis reconstruction 2-D... ")
#
# for l in 1:5
#
#     pos_coeffs = standard_coeffs(x->sin(4*x[1]+x[2]),(l,l))
#     pos_err = x->(standard_reconstruct(pos_coeffs, (l,l), (x[1],x[2]))-sin(4*x[1]+x[2]))^2
#     @test hcubature(pos_err, [0,0], [1,1]; abstol=1.0e-9, maxevals=500)[1]<1/(1<<(l))
# end
#
# println("Test Passed.")
#
# print("Testing full hat basis reconstruction 2-D... ")
#
# for l in 3:7
#     hier_coeffs = coeffs_hat(2, l, x->sin(4*x[1]+x[2]); scheme="full")
#     hier_err = x->(reconstruct_hat(hier_coeffs, (x[1],x[2]))-sin(4*x[1]+x[2]))^2
#     @test hcubature(hier_err, [0,0], [1,1]; abstol=1.0e-9)[1]<1/(1<<(l))
# end
#
# println("Test Passed.")
#
# print("Testing sparse hat basis reconstruction 2-D... ")
#
# for l in 1:5
#     sparse_coeffs = coeffs_hat(2, l, x->sin(4*x[1]+x[2]); scheme="sparse")
#     sparse_err = x->(reconstruct_hat(sparse_coeffs, (x[1],x[2]))-sin(4*x[1]+x[2]))^2
#     @test hcubature(sparse_err, [0,0], [1,1];abstol=1.0e-9)[1]<1/(1<<(l))
# end
#
# println("Test Passed.")
#
# #---------------------------------------
# # Testing regular hier DG reconstruction
# #---------------------------------------
#
# print("Testing full discontinuous Galerkin (DG) reconstruction 1-D... ")
#
# for k in 1:5
#     for l in 1:6
#         dict = coeffs_DG(1, k, l, x->sin(4*x[1]); scheme="full")
#         err = x->(reconstruct_DG(dict, [x])-sin(4*x))^2
#         @test hquadrature(err, 0, 1; abstol=1.0e-10)[1]< 1/(1<<(l+k-1))
#     end
# end
#
# for k in 1:5
#     for l in 1:6
#         dict = coeffs_DG(1, k, l, x->sin(4*x[1]); scheme="sparse")
#         err = x->(reconstruct_DG(dict, [x])-sin(4*x))^2
#         @test hquadrature(err, 0, 1; abstol=1.0e-10)[1]< 1/(1<<(l+k-1))
#     end
# end
#
# println("Test Passed.")
#
# # Sparsification should be trivial in 1D
#
# print("Testing that both reconstructions are equivalent... ")
#
# for k in 1:5
#     for l in 1:6
#         full_coeffs = coeffs_DG(1, k, l, x->sin(4*x[1]); scheme="full")
#         sparse_coeffs = coeffs_DG(1, k, l, x->sin(4*x[1]); scheme="full")
#         diff_full = x->(reconstruct_DG(full_coeffs, [x])-sin(4*x))^2
#         diff_sparse = x->(reconstruct_DG(sparse_coeffs, [x])-sin(4*x))^2
#         diff_both = x->(reconstruct_DG(full_coeffs, [x])-reconstruct_DG(sparse_coeffs, [x]))^2
#
#         @test abs(hquadrature(diff_full,0,1;abstol=1.0e-9)[1]-
#                   hquadrature(diff_sparse,0,1;abstol=1.0e-9)[1])<1.0e-14
#         @test hquadrature(diff_both,0,1;abstol=1.0e-9)[1]<1.0e-14
#     end
# end
#
# println("Test Passed.")
#
# # Testing sparse DG reconstruction in multidimensional space:
#
# print("Testing sparse DG basis reconstruction 2-D... ")
#
# for k in 1:5
#     for l in 1:6
#         dict = coeffs_DG(2, k, l, x->sin(4*x[1]+x[2]))
#         err = x->(reconstruct_DG(dict, [x[1],x[2]])-sin(4*x[1]+x[2]))^2
#         @test hcubature(err, [0,0], [1,1]; abstol=1.0e-10, maxevals=500)[1]< 1/(1<<(l+k-2))
#     end
# end
#
# println("Test Passed.")
#
# #---------------------------------------
# # Testing vector hier DG reconstruction
# #---------------------------------------
#
# print("Testing full DG reconstruction 1-D with vector coefficients... ")
#
# for k in 1:5
#     for l in 1:6
#         vect = vcoeffs_DG(1, k, l, x->sin(4*x[1]); scheme="full")
#         dict = V2D(1, k, l, vect; scheme="full")
#         err = x->(reconstruct_DG(dict,[x[1]])-sin(4*x[1]))^2
#         @test hquadrature(err,0,1; reltol=1.0e-9, abstol=1.0e-12)[1] < 1/(1<<((l+k-1)))
#     end
# end
#
# println("Test Passed.")
#
# # Testing D2V and V2D
#
# # Full case:
#
# print("Testing D2V and V2D Full Case 2D... ")
#
# for k in 1:5
#     for l in 1:4
#         dict = coeffs_DG(2, k, l, x->sin(4*x[1]+x[2]); scheme="full")
#         vect = D2V(2, k, l, dict; scheme="full")
#         @test V2D(2, k, l, vect; scheme="full")==dict
#     end
# end
#
# for k in 1:5
#     for l in 1:4
#         vect = vcoeffs_DG(2, k, l, x->sin(4*x[1]+x[2]); scheme="full")
#         dict = V2D(2, k, l, vect; scheme="full")
#         @test D2V(2, k, l, dict; scheme="full")==vect
#     end
# end
#
# println("Test Passed.")
#
# # Sparse case:
#
# print("Testing D2V and V2D Sparse Case 2D... ")
#
# for k in 1:5
#     for l in 1:4
#         dict = coeffs_DG(2, k, l, x->sin(4*x[1]+x[2]); scheme="sparse")
#         vect = D2V(2, k, l, dict; scheme="sparse")
#         @test V2D(2, k, l, vect; scheme="sparse")==dict
#     end
# end
#
# for k in 1:5
#     for l in 1:4
#         vect = vcoeffs_DG(2, k, l, x->sin(4*x[1]+x[2]); scheme="sparse")
#         dict= V2D(2, k, l, vect; scheme="sparse")
#         @test D2V(2, k, l, dict; scheme="sparse")==vect
#     end
# end
#
# println("Test Passed.")
#
# #--------------------------------------
# # Testing Differentiation
# #--------------------------------------
#
# print("Testing differentiation 1-D DG basis... ")
#
# k=3
# for l in 2:5
#     frefVD = V2Dref(1, k, l);
#     frefDV = D2Vref(1, k, l);
#     D_op = D_matrix(1, k, l, frefVD, frefDV; scheme="full")
#     vcoeffs = vcoeffs_DG(1, k, l, x->cos(2*pi*x[1]); scheme="full")
#     dvcoeffs = *(D_op,vcoeffs)
#     dict= V2D(1, k, l, vcoeffs; scheme="full")
#     ddict = V2D(1, k, l, dvcoeffs; scheme="full")
#     err = x->(reconstruct_DG(ddict,[x[1]])+2*pi*sin(2*pi*x[1]))^2
#     @test hquadrature(err, 0, 1, abstol=1.0e-10, maxevals=500)[1]<1/(1<<(k+l-2))
# end
#
# println("Test Passed.")
#
# print("Testing differentiation 2-D full DG basis... ")
#
# D = 2
# k = 3
# for l in 2:5
#     D_op = D_matrix(D, 1, k, l; scheme="full")
#     vcoeffs = vcoeffs_DG(D, k, l, x->cos(2*pi*x[1])*cos(2*pi*x[2]); scheme="full")
#     dvcoeffs = *(D_op, vcoeffs)
#     dict= V2D(D, k, l, vcoeffs; scheme="full")
#     ddict = V2D(D, k, l, dvcoeffs; scheme="full")
#     err = x->(reconstruct_DG(ddict,[x[1],x[2]])+2*pi*sin(2*pi*x[1])*cos(2*pi*x[2]))^2
#     @test hcubature(err, [0,0], [1,1], abstol=1.0e-10, maxevals=500)[1] < 1/(1<<(k+l-2))
# end
#
# println("Test Passed.")
#

#--------------------------------------
# Testing PDE solvers
#--------------------------------------

# The wave equation, testing conservation of energy

# In the position basis:

import GalerkinSparseGrids.energy_func_1D

print("Testing wave equation solver 1-D position DG basis... ")

k = 4
level = 4
f0 = x->sin(2*pi*x[1])
v0 = x->2*pi*cos(2*pi*x[1])
pos_soln = wave_evolve_1D(k, level, f0, v0, 0, 1; base="pos", order="45")
energy_soln = energy_func_1D(k, level, pos_soln; base="pos")
for i in 1:length(energy_soln[2])
    @test abs(sqrt(energy_soln[2][i])-2*pi)<=1.0e-7
end

println("Test Passed.")

# In the hierarchical basis:

print("Testing wave equation solver 1-D hierarchical DG basis... ")

hier_soln  = wave_evolve_1D(k, level, f0, v0, 0, 1; base="hier", order="45")
energy_soln = energy_func_1D(k, level, hier_soln; base="hier")
for i in 1:length(energy_soln[2])
    @test abs(sqrt(energy_soln[2][i])-2*pi)<=1.0e-7
end

println("Test Passed.")

# 2D Sparse case:

import GalerkinSparseGrids.energy_func

print("Testing wave equation solver on 2-D sparse DG basis... ")

D = 2
k_used = 3
n_used = 5
f0 = x->sin(2*pi*x[1])*sin(2*pi*x[2])
v0 = x->0
sparse_soln = wave_evolve(D, k_used, n_used, f0, v0, 0,1; order = "78", scheme = "sparse");
senergy = energy_func(D, k_used, n_used, sparse_soln)
@test senergy[2][1]-senergy[2][end]>0
@test abs(senergy[2][1]-senergy[2][end])<1.0e-8
for i in 1:length(senergy[1])
    @test abs(sqrt(senergy[2][i])-sqrt(2)*pi)<=1.0e-4
end

sparse_soln = wave_evolve(D, k_used, n_used, f0, v0, 0,1; order = "45", scheme = "sparse");
senergy = energy_func(D, k_used, n_used, sparse_soln)
@test senergy[2][1]-senergy[2][end]>0
@test abs(senergy[2][1]-senergy[2][end])<1.0e-8
for i in 1:length(senergy[1])
    @test abs(sqrt(senergy[2][i])-sqrt(2)*pi)<=1.0e-4
end

println("Test Passed.")

print("Testing traveling wave example in 2-D... ")

m = [1,2]
truesoln = x -> cos(2*pi*(vecdot(m,x) - sqrt(vecdot(m,m))*0.54))

D = length(m)
k_used = 3
n_used = 6

soln = traveling_wave_solver(k_used, n_used, m, 0, 0.54; order="45")
dict = V2D(D, k_used, n_used, soln[2][end]; scheme="sparse")
@test mcerr(x->reconstruct_DG(dict, [x...]), truesoln, D) < 0.05

soln = traveling_wave_solver(k_used, n_used, m, 0, 0.54; order="78")
dict = V2D(D, k_used, n_used, soln[2][end]; scheme="sparse")
@test mcerr(x->reconstruct_DG(dict, [x...]), truesoln, D) < 0.05

println("Test passed.")

println("Finished Tests of GalerkinSparseGrids.jl---------------------------------------")
