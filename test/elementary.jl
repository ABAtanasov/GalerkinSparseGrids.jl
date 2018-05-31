using GalerkinSparseGrids
using Base.Test
using Cubature

#--------------------------------------
# Elementary Tests
#--------------------------------------

ε = eps(Float64)

print("Testing cell_index function... ")
import GalerkinSparseGrids.cell_index
for l in 1:5
	@test cell_index(1.3,l)==(1<<(l-1))
end
for l in 1:5
	@test cell_index(.01,l)==1
end
@test cell_index(.3,1)==1
@test cell_index(.3,2)==1
@test cell_index(.3,3)==2
@test cell_index(.3,4)==3
@test cell_index(.3,5)==5
println("Tests Passed.")

import GalerkinSparseGrids.inner_product

print("Testing 1D inner product... ")
@test inner_product(x->x[1]^2, x->x[1]^3, (0,),CartesianIndex((1,))) ≈ (1/6) atol=2*ε
@test inner_product(x->sin(pi*x[1]), x->cos(pi*x[1]), (0,),CartesianIndex((1,))) ≈ 0 atol=2*ε
println("Tests Passed.")

print("Testing 2D inner product... ")
@test inner_product(x->(x[1]^2+x[2]^2), x->x[1]^3, (0,0),CartesianIndex((1,1))) ≈ .25 atol=2*ε
println("Tests Passed.")

print("Testing 3D inner product... ")
@test inner_product(x->(x[1]^2+x[2]^2-x[3]^2), x->x[1]^3+x[2]+x[3], (0,0,0),CartesianIndex((1,1,1))) ≈ .5 atol=2*ε
println("Tests Passed.")
