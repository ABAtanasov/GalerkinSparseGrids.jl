using GalerkinSparseGrids
using Base.Test
using Cubature
using ODE


println("Beginning Tests of GalerkinSparseGrids.jl--------------------------------------")

#--------------------------------------
# Elementary Tests
#--------------------------------------

import GalerkinSparseGrids.pos
print("Testing pos function... ")
for i in -10:0
	@test pos(i)==0
end
for i in 0:10
	@test pos(i)==i
end
println("Tests Passed.")

print("Testing hat_index function... ")
import GalerkinSparseGrids.hat_index
for l in 1:5
	@test hat_index(1.3,l)==(1<<(l-1))
end
for l in 1:5
	@test hat_index(.01,l)==1
end
@test hat_index(.3,1)==1
@test hat_index(.3,2)==1
@test hat_index(.3,3)==2
@test hat_index(.3,4)==3
@test hat_index(.3,5)==5
println("Tests Passed.")

import GalerkinSparseGrids.inner_product
print("Testing 1D inner product... ")
@test abs(inner_product(x->x[1]^2, x->x[1]^3, (0,),CartesianIndex((1,)))-(1/6)) < 2*eps(Float64)
@test abs(inner_product(x->sin(pi*x[1]), x->cos(pi*x[1]), (0,),CartesianIndex((1,)))) < eps(Float64)
println("Tests Passed.")
print("Testing 2D inner product... ")
@test abs(inner_product(x->(x[1]^2+x[2]^2), x->x[1]^3, (0,0),CartesianIndex((1,1)))-.25)<2*eps(Float64)
println("Tests Passed.")
print("Testing 3D inner product... ")
@test abs(inner_product(x->(x[1]^2+x[2]^2-x[3]^2), x->x[1]^3+x[2]+x[3], (0,0,0),CartesianIndex((1,1,1)))-.5)<2*eps(Float64)
println("Tests Passed.")

#--------------------------------------
# Testing hat reconstruction
#--------------------------------------

# Position:
print("Testing position hat basis reconstruction 1-D... ")

for l in 1:7
    pos_coeffs = standard_coefficients(x->sin(4*x[1]),(l,))
    @test hquadrature(x->(standard_reconstruct(pos_coeffs, (l,), (x,))-sin(4*x))^2,0,1;abstol=1.0e-9)[1]<1/(1<<(l))
end

println("Test Passed.")

# Hierarchical:

print("Testing hierarchical hat basis reconstruction 1-D... ")

for l in 3:9
    hier_coeffs = hier_coefficients(x->sin(4*x[1]),(l,))
    @test hquadrature(x->(reconstruct(hier_coeffs, (x,))-sin(4*x))^2,0,1;abstol=1.0e-9)[1]<1/(1<<(l))
end

println("Test Passed.")

# Hierarchical and position basis should yield the same results for hat functions:

print("Testing that both reconstructions are equivalent... ")

for l in 1:7
    pos_coeffs=standard_coefficients(x->sin(4*x[1]),(l,))
    hier_coeffs = hier_coefficients(x->sin(4*x[1]),(l+2,))
    @test abs(hquadrature(x->(standard_reconstruct(pos_coeffs, (l,), (x,))-sin(4*x))^2,0,1;abstol=1.0e-9)[1]-
            hquadrature(x->(reconstruct(hier_coeffs, (x,))-sin(4*x))^2,0,1;abstol=1.0e-9)[1])<1.0e-14
end

println("Test Passed.")

# Multidimensional hat reconstruction:

print("Testing position hat basis reconstruction 2-D... ")

for l in 1:5

    hier_coeffs = standard_coefficients(x->sin(4*x[1]+x[2]),(l,l))
    @test hcubature(x->(standard_reconstruct(hier_coeffs, (l,l), (x[1],x[2]))-sin(4*x[1]+x[2]))^2,[0,0],[1,1];abstol=1.0e-9, maxevals=500)[1]<1/(1<<(l))
end

println("Test Passed.")

print("Testing hierarchical hat basis reconstruction 2-D... ")

for l in 3:7
    hier_coeffs = hier_coefficients(x->sin(4*x[1]+x[2]),(l,l))
    @test hcubature(x->(reconstruct(hier_coeffs, (x[1],x[2]))-sin(4*x[1]+x[2]))^2,[0,0],[1,1];abstol=1.0e-9)[1]<1/(1<<(l))
end

println("Test Passed.")

print("Testing sparse hat basis reconstruction 2-D... ")

for l in 1:5
    sparse_coeffs = sparse_coefficients(x->sin(4*x[1]+x[2]),l,2)
    @test hcubature(x->(reconstruct(sparse_coeffs, (x[1],x[2]))-sin(4*x[1]+x[2]))^2,[0,0],[1,1];abstol=1.0e-9)[1]<1/(1<<(l))
end

println("Test Passed.")


#---------------------------------------
# Testing regular hier DG reconstruction
#---------------------------------------

print("Testing hierarchical discontinuous Galerkin (DG) basis reconstruction 1-D... ")

for k in 1:5
    for l in 1:6
        dict = hier_coefficients_DG(k,x->sin(4*x[1]),(l,))
        @test hquadrature(x->(reconstruct_DG(k, dict, [x])-sin(4*x))^2, 0, 1; abstol=1.0e-10)[1]< 1/(1<<(l+k-1))
    end
end

for k in 1:5
    for l in 1:6
        dict = sparse_coefficients_DG(k,x->sin(4*x[1]),l,1)
        @test hquadrature(x->(reconstruct_DG(k, dict, [x])-sin(4*x))^2, 0, 1; abstol=1.0e-10)[1]< 1/(1<<(l+k-1))
    end
end

println("Test Passed.")

# Testing sparse DG reconstruction in multidimensional space:

print("Testing sparse DG basis reconstruction 2-D... ")

for k in 1:5
    for l in 1:6
        dict = sparse_coefficients_DG(k,x->sin(4*x[1]+x[2]),l,2)
        @test hcubature(x->(reconstruct_DG(k, dict, [x[1],x[2]])-sin(4*x[1]+x[2]))^2, [0,0], [1,1]; abstol=1.0e-10, maxevals=500)[1]< 1/(1<<(l+k-2))
    end
end

println("Test Passed.")

# Testing vector reconstruction:

print("Testing hierarchical DG basis reconstruction 1-D with vector coefficients... ")

for k in 1:5
    for l in 1:6
        vect = vhier_coefficients_DG(k, x->sin(4*x[1]), (l,))
        dict = full_V2D(k,vect,(l,))
        @test hquadrature(x-> (reconstruct_DG(k,dict,[x[1]])-sin(4*x[1]))^2,0,1; reltol=1.0e-9, abstol=1.0e-12)[1] < 1/(1<<((l+k-1)))
    end
end

println("Test Passed.")

# Testing D2V and V2D

# Full case:

print("Testing D2V and V2D Full Case 2D... ")

for k in 1:5
    for l in 1:5
        dict = hier_coefficients_DG(k,x->sin(4*x[1]+x[2]),(l,l))
        vect = full_D2V(k,dict,(l,l))
        @test full_V2D(k,vect,(l,l))==dict
    end
end

for k in 1:5
    for l in 1:5
        vect = vhier_coefficients_DG(k,x->sin(4*x[1]+x[2]),(l,l))
        dict = full_V2D(k,vect,(l,l))
        @test full_D2V(k,dict,(l,l))==vect
    end
end

println("Test Passed.")

# Sparse case:

print("Testing D2V and V2D Sparse Case 2D... ")

for k in 1:5
    for l in 1:5
        dict = sparse_coefficients_DG(k,x->sin(4*x[1]+x[2]),l,2)
        vect = sparse_D2V(k,dict,l)
        @test sparse_V2D(k,vect,l,2)==dict
    end
end

for k in 1:5
    for l in 1:5
        vect = vsparse_coefficients_DG(k,x->sin(4*x[1]+x[2]),l,2)
        dict= sparse_V2D(k,vect,l,2)
        @test sparse_D2V(k,dict,l)==vect
    end
end

println("Test Passed.")

#--------------------------------------
# Testing Differentiation 
#--------------------------------------

print("Testing differentiation 1-D DG basis... ")

k=3
for l in 2:5
	levels = (l+1,)
    frefVD = full_referenceV2D(k, levels);
    frefDV = full_referenceD2V(k, levels);
    D_op = full_D_matrix(1, k, l, frefVD, frefDV)
    vcoeffs = vhier_coefficients_DG(k, x->cos(2*pi*x[1]), levels)
    dvcoeffs = *(D_op,vcoeffs)
    dict= full_V2D(k, vcoeffs, levels)
    ddict = full_V2D(k, dvcoeffs, levels)
    err = hquadrature(x->(reconstruct_DG(k,ddict,[x[1]])+2*pi*sin(2*pi*x[1]))^2, 0, 1, abstol=1.0e-10, maxevals=500)[1]
    @test err<1/(1<<(k+l-2))
end

println("Test Passed.")

print("Testing differentiation 2-D full DG basis... ")

k=3
for l in 2:5
	levels = (l+1, l+1)
    frefVD = full_referenceV2D(k, levels);
    frefDV = full_referenceD2V(k, levels);
    D_op = full_D_matrix(1, k, l, frefVD, frefDV)
    vcoeffs = vhier_coefficients_DG(k, x->cos(2*pi*x[1])*cos(2*pi*x[2]), levels)
    dvcoeffs = *(D_op,vcoeffs)
    dict= full_V2D(k,vcoeffs, levels)
    ddict = full_V2D(k,dvcoeffs, levels)
    err = hcubature(x->(reconstruct_DG(k,ddict,[x[1],x[2]])+2*pi*sin(2*pi*x[1])*cos(2*pi*x[2]))^2,[0,0],[1,1],abstol=1.0e-10,maxevals=500)[1]
    @test err<1/(1<<(k+l-2))
end

println("Test Passed.")

#--------------------------------------
# Testing PDE solvers
#--------------------------------------

# The wave equation, testing conservation of energy

# In the position basis:

import GalerkinSparseGrids.pos_energy_func

print("Testing wave equation solver 1-D position DG basis... ")

pos_soln=pos_wave_equation45(x->sin(2*pi*x),x->2*pi*cos(2*pi*x), 4,4,0,1);
energy_soln=pos_energy_func(4,4,pos_soln)
for i in 1:length(energy_soln[2])
    @test abs(sqrt(energy_soln[2][i])-2*pi)<=1.0e-7
end

println("Test Passed.")

# In the hierarchical basis:

print("Testing wave equation solver 1-D hierarchical DG basis... ")

import GalerkinSparseGrids.hier_energy_func

hier_soln=hier_wave_equation45(x->sin(2*pi*x[1]),x->2*pi*cos(2*pi*x[1]), 4,4,0,1);
energy_soln=hier_energy_func(4,4,hier_soln)
for i in 1:length(energy_soln[2])
    @test abs(sqrt(energy_soln[2][i])-2*pi)<=1.0e-7
end

println("Test Passed.")


import GalerkinSparseGrids.sparse_energy_func

print("Testing wave equation solver on 2-D sparse DG basis... ")

sparse_soln=sparse_wave_equation78(x->sin(2*pi*x[1])*sin(2*pi*x[2]),x->0,3,5,2,0,1);
senergy=sparse_energy_func(3,5,2,sparse_soln)
@test senergy[2][1]-senergy[2][end]>0
@test abs(senergy[2][1]-senergy[2][end])<1.0e-8

for i in 1:length(senergy[1])
    @test abs(sqrt(senergy[2][i])-sqrt(2)*pi)<=1.0e-4
end

soln=sparse_wave_equation45(x->sin(2*pi*x[1])*sin(2*pi*x[2]),x->0,3,5,2,0,1);
senergy=sparse_energy_func(3,5,2,soln)
@test senergy[2][1]-senergy[2][end]>0
@test abs(senergy[2][1]-senergy[2][end])<1.0e-8

for i in 1:length(senergy[1])
    @test abs(sqrt(senergy[2][i])-sqrt(2)*pi)<=1.0e-4
end

println("Test Passed.")

print("Testing traveling wave example in 2-D... ")

m = [1,2]
truesoln = x -> cos(2*pi*(vecdot(m,x) - sqrt(vecdot(m,m))*0.54))

k_used = 3
n_used = 6

D = length(m)

soln = traveling_wave_equation45(k_used, n_used, m, 0, 0.54)
dict = sparse_V2D(k_used, soln[2][end], n_used, D)

@test mcerr(x->reconstruct_DG(k_used, dict, [x...]), truesoln, D) < 0.05

println("Test passed.")

println("Finished Tests of GalerkinSparseGrids.jl---------------------------------------")
