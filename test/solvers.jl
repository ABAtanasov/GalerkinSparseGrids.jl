using GalerkinSparseGrids
using Base.Test
using Cubature
using ODE

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
sparse_soln = wave_evolve(D, k_used, n_used, f0, v0, 0,1; order="78", scheme="sparse");
senergy = energy_func(D, k_used, n_used, sparse_soln)
@test senergy[2][1]-senergy[2][end] > 0
@test abs(senergy[2][1]-senergy[2][end]) < 1.0e-8
for i in 1:length(senergy[1])
    @test abs(sqrt(senergy[2][i])-sqrt(2)*pi) <= 1.0e-4
end

sparse_soln = wave_evolve(D, k_used, n_used, f0, v0, 0,1; order="45", scheme="sparse");
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
