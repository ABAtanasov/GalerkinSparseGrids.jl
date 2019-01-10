using GalerkinSparseGrids
using Test
using HCubature
using StaticArrays
using ODE
using SparseArrays

#--------------------------------------
# Testing PDE solvers
#--------------------------------------

# The wave equation, testing conservation of energy

# In the position basis:

import GalerkinSparseGrids.energy_func_1D
import GalerkinSparseGrids.energy_func

@testset "solvers.jl" begin
    @info "Testing wave equation solver 1-D position DG basis... "

    global k = 4
    global level = 4
    global f0 = x->sin(2*pi*x[1])
    global v0 = x->2*pi*cos(2*pi*x[1])
    pos_soln = wave_evolve_1D(k, level, f0, v0, 0, 1; basis="pos", order="45") #problem is here
    energy_soln = energy_func_1D(k, level, pos_soln; basis="pos")
    @testset "Wave equation solver 1D position DG basis" begin
        for i in 1:length(energy_soln[2])
            @test sqrt(energy_soln[2][i]) ≈ 2*pi atol=1.0e-7
        end
    end

    # In the hierarchical basis:
    @info "Testing wave equation solver 1-D hierarchical DG basis... "
    @testset "Wave equation solver 1D hier DG basis" begin
        hier_soln  = wave_evolve_1D(k, level, f0, v0, 0, 1; basis="hier", order="45")
        energy_soln = energy_func_1D(k, level, hier_soln; basis="hier")
        for i in 1:length(energy_soln[2])
            @test sqrt(energy_soln[2][i]) ≈ 2*pi atol=1.0e-7
        end

        hier_soln  = wave_evolve(1, k, level, f0, v0, 0, 1; order="78")
        energy_soln = energy_func(1, k, level, hier_soln)
        for i in 1:length(energy_soln[2])
            @test sqrt(energy_soln[2][i]) ≈ 2*pi atol=1.0e-7
        end
    end

    # 2D Sparse case:

    @info "Testing wave equation solver on 2-D sparse DG basis... "

    @testset "Wave equation solver 2D sparse DG basis" begin

        D = 2
        k_used = 3
        n_used = 5

        f0 = x->sin(2*pi*x[1])*sin(2*pi*x[2])
        v0 = x->0
        sparse_soln = wave_evolve(D, k_used, n_used, f0, v0, 0, 1; order="78", scheme="sparse");
        senergy = energy_func(D, k_used, n_used, sparse_soln)
        @test senergy[2][1]-senergy[2][end] > 0
        @test abs(senergy[2][1]-senergy[2][end]) < 1.0e-8
        for i in 1:length(senergy[1])
            @test sqrt(senergy[2][i]) ≈ sqrt(2)*pi atol=1.0e-4
        end

        sparse_soln = wave_evolve(D, k_used, n_used, f0, v0, 0,1; order="45", scheme="sparse");
        senergy = energy_func(D, k_used, n_used, sparse_soln)
        @test senergy[2][1] > senergy[2][end]
        @test senergy[2][1] ≈ senergy[2][end] atol=1.0e-8
        for i in 1:length(senergy[1])
            @test sqrt(senergy[2][i]) ≈ sqrt(2)*pi atol=1.0e-4
        end
    end

    # @info "Testing traveling wave example in 2-D... "
    # @testset "Traveling wave example 2D" begin
    #     m = [1, 2]
    #     truesoln = x -> cos(2*pi*(dot(m,x) - sqrt(dot(m,m))*0.54))
    #
    #     D = 2
    #     k_used = 3
    #     n_used = 6
    #
    #     soln = traveling_wave_solver(k_used, n_used, m, 0, 0.54; order="45")
    #     dict = V2D(D, k_used, n_used, soln[2][end]; scheme="sparse")
    #     @test mcerr(x->reconstruct_DG(dict, [x...]), truesoln, D) < 0.05
    #
    #     soln = traveling_wave_solver(k_used, n_used, m, 0, 0.54; order="78")
    #     dict = V2D(D, k_used, n_used, soln[2][end]; scheme="sparse")
    #     @test mcerr(x->reconstruct_DG(dict, [x...]), truesoln, D) < 0.05
    # end
end
