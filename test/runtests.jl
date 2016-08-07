using GalerkinSparseGrids
using Base.Test
using Cubature
using ODE


#--------------------------------------
# Elementary Tests
#--------------------------------------

import GalerkinSparseGrids.pos
for i in -10:0
	@test pos(i)==0
end
for i in 0:10
	@test pos(i)==i
end

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


import GalerkinSparseGrids.inner_product
print("Testing 1D inner product... ")
@test abs(inner_product(x->x[1]^2, x->x[1]^3, (0,),CartesianIndex((1,)))-(1/6)) < 2*eps(Float64)
@test abs(inner_product(x->sin(pi*x[1]), x->cos(pi*x[1]), (0,),CartesianIndex((1,)))) < eps(Float64)
println("Test Passed.")
print("Testing 2D inner product... ")
@test abs(inner_product(x->(x[1]^2+x[2]^2), x->x[1]^3, (0,0),CartesianIndex((1,1)))-.25)<2*eps(Float64)
println("Test Passed.")
print("Testing 3D inner product... ")
@test abs(inner_product(x->(x[1]^2+x[2]^2-x[3]^2), x->x[1]^3+x[2]+x[3], (0,0,0),CartesianIndex((1,1,1)))-.5)<2*eps(Float64)
println("Test Passed.")

#--------------------------------------
# Testing hat reconstruction 
#--------------------------------------

#position:
print("Testing position hat basis reconstruction 1-D... ")

for l in 1:7
    @test hquadrature(x->(standard_reconstruct(standard_coefficients(x->sin(4*x[1]),(l,)), (l,), (x,))-sin(4*x))^2,0,1;abstol=1.0e-9)[1]<1/(1<<(l))
end

println("Test Passed.")

#hierarchical:

print("Testing hierarchical hat basis reconstruction 1-D... ")

for l in 3:9
    @test hquadrature(x->(reconstruct(hier_coefficients(x->sin(4*x[1]),(l,)), (x,))-sin(4*x))^2,0,1;abstol=1.0e-9)[1]<1/(1<<(l))
end

println("Test Passed.")

#hierarchical and position basis should yield the same results for hat functions:

print("Testing that both reconstructions are equivalent... ")

for l in 1:7
    @test abs(hquadrature(x->(standard_reconstruct(standard_coefficients(x->sin(4*x[1]),(l,)), (l,), (x,))-sin(4*x))^2,0,1;abstol=1.0e-9)[1]-
			hquadrature(x->(reconstruct(hier_coefficients(x->sin(4*x[1]),(l+2,)), (x,))-sin(4*x))^2,0,1;abstol=1.0e-9)[1])<1.0e-14
end

println("Test Passed.")



#multidimensional hat reconstruction: 


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

# testing sparse DG reconstruction in multidimensional space

print("Testing sparse DG basis reconstruction 2-D... ")

for k in 1:5
    for l in 1:6
        dict = sparse_coefficients_DG(k,x->sin(4*x[1]+x[2]),l,2)
        @test hcubature(x->(reconstruct_DG(k, dict, [x[1],x[2]])-sin(4*x[1]+x[2]))^2, [0,0], [1,1]; abstol=1.0e-10, maxevals=500)[1]< 1/(1<<(l+k-2))     
    end
end

println("Test Passed.")

# testing vector reconstruction 

print("Testing hierarchical DG basis reconstruction 1-D with vector coefficients... ")

for k in 1:5
    for l in 1:6
        vect = vhier_coefficients_DG(k, x->sin(4*x[1]), (l,))
        dict = Full_V2D(k,vect,(l,))
        @test hquadrature(x-> (reconstruct_DG(k,dict,[x[1]])-sin(4*x[1]))^2,0,1; reltol=1.0e-9, abstol=1.0e-12)[1] < 1/(1<<((l+k-1)))
    end
end

println("Test Passed.")

#--------------------------------------
# testing differentiation 
#--------------------------------------

print("Testing differentiation 1-D DG basis... ")

for k in 1:5
    for l in 1:6
	    vect = vhier_coefficients_DG(3, x->sin(4*x[1]), (5,))
	    refVD=full_referenceV2D(3,(5,)); 
	    refDV=full_referenceD2V(3,(5,));
	    sD = sD_matrix(1,3, refVD, refDV)
	    dvect = *(sD, vect)
	    dict=Full_V2D(3,dvect,(5,))
	    @test hquadrature(x->(reconstruct_DG(3,dict,[x[1]])-4*cos(4*x[1]))^2,0,1; abstol=1.0e-10)[1] < 1/(1<<((l+k-2)))
    end
end

println("Test Passed.")

#--------------------------------------
# testing PDE solvers
#--------------------------------------

#a wave equation, testing conservation of energy

#in position basis
print("Testing wave equation solver 1-D position DG basis... ")

pos_soln=pos_wave_equation45(x->sin(2*pi*x),x->2*pi*cos(2*pi*x), 4,4,0,1);
energy_soln=pos_energy_func(4,4,pos_soln)
for i in 1:length(energy_soln[2])
    @test abs(sqrt(energy_soln[2][i])-2*pi)<=1.0e-7
end

println("Test Passed.")

#in hierarchical basis
print("Testing wave equation solver 1-D hierarchical DG basis... ")

pos_soln=hier_wave_equation45(x->sin(2*pi*x[1]),x->2*pi*cos(2*pi*x[1]), 4,4,0,1);
energy_soln=hier_energy_func(4,4,pos_soln)
for i in 1:length(energy_soln[2])
    @test abs(sqrt(energy_soln[2][i])-2*pi)<=1.0e-7
end

println("Test Passed.")
