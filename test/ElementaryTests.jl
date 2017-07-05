using GalerkinSparseGrids
using Base.Test
using Cubature

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
