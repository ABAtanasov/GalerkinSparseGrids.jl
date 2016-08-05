using GalerkinSparseGrids
using Base.Test
using Cubature
using ODE

#--------------------------------------
# Testing hat reconstruction 
#--------------------------------------

#position:
for l in 1:7
    @test hquadrature(x->(standard_reconstruct(standard_coefficients(x->sin(4*x[1]),(l,)), (l,), (x,))-sin(4*x))^2,0,1;abstol=1.0e-9)[1]<1/(1<<(l))
end

#hierarchical:
for l in 3:9
    @test hquadrature(x->(reconstruct(hier_coefficients(x->sin(4*x[1]),(l,)), (x,))-sin(4*x))^2,0,1;abstol=1.0e-9)[1]<1/(1<<(l))
end


#hierarchical and position basis should yield the same results for hat functions:
for l in 1:7
    @test abs(hquadrature(x->(standard_reconstruct(standard_coefficients(x->sin(4*x[1]),(l,)), (l,), (x,))-sin(4*x))^2,0,1;abstol=1.0e-9)[1]-
			hquadrature(x->(reconstruct(hier_coefficients(x->sin(4*x[1]),(l+2,)), (x,))-sin(4*x))^2,0,1;abstol=1.0e-9)[1])<1.0e-14
end

#multidimensional hat reconstruction: 


#---------------------------------------
# Testing regular hier DG reconstruction 
#---------------------------------------

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

# testing sparse DG reconstruction in multidimensional space

for k in 1:5
    for l in 1:6
        dict = sparse_coefficients_DG(k,x->sin(4*x[1]+x[2]),l,2)
        @test hcubature(x->(reconstruct_DG(k, dict, [x[1],x[2]])-sin(4*x[1]+x[2]))^2, [0,0], [1,1]; abstol=1.0e-10, maxevals=500)[1]< 1/(1<<(l+k-2))     
    end
end

# testing vector reconstruction 

for k in 1:5
    for l in 1:6
        vect = vhier_coefficients_DG(k, x->sin(4*x[1]), (l,))
        dict = Full_V2D(k,vect,(l,))
        @test hquadrature(x-> (reconstruct_DG(k,dict,[x[1]])-sin(4*x[1]))^2,0,1; reltol=1.0e-9, abstol=1.0e-12)[1] < 1/(1<<((l+k-1)))
    end
end

#--------------------------------------
# testing differentiation 
#--------------------------------------

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

#--------------------------------------
# testing PDE solvers
#--------------------------------------

#a wave equation, testing conservation of energy

#in position basis
pos_soln=pos_wave_equation45(x->sin(2*pi*x),x->2*pi*cos(2*pi*x), 4,4,0,1);
energy_soln=pos_energy_func(4,4,pos_soln)
for i in 1:length(energy_soln[2])
    @test abs(sqrt(energy_soln[2][i])-2*pi)<=1.0e-7
end

#in hierarchical basis
pos_soln=hier_wave_equation45(x->sin(2*pi*x[1]),x->2*pi*cos(2*pi*x[1]), 4,4,0,1);
energy_soln=hier_energy_func(4,4,pos_soln)
for i in 1:length(energy_soln[2])
    @test abs(sqrt(energy_soln[2][i])-2*pi)<=1.0e-7
end

