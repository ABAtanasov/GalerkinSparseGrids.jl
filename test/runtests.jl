using GalerkinSparseGrids
using Base.Test
using Cubature

@test 1 == 1

for k in 1:5
    for l in 1:6
        vect = vhier_coefficients_DG(k, x->sin(4*x[1]), (l,))
        dict = Full_V2D(k,vect,(l,))
        @test hquadrature(x-> (reconstruct_DG(k,dict,[x[1]])-sin(4*x[1]))^2,0,1; 
        reltol=1.0e-9, abstol=1.0e-12)[1] < 1/(1<<((l+k-1)))
    end
end


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

