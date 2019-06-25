using GalerkinSparseGrids
using Test
using HCubature
using ODE
using SparseArrays
using LinearAlgebra

tests = [
    "elementary.jl",
    "basis_constructor.jl",
    "hier_DG.jl",
    "vhier_DG.jl",
    "transformations.jl",
    "differentiation.jl",
    "solvers.jl"
]

@testset "GalerkinSparseGrids" begin
    for filename in tests
        name = first(splitext(filename))
        include(filename)
    end
end

# Later, we will also test the examples:

# include(joinpath(dirname(@__FILE__), "..", "examples", "interpolation.jl"))
# include(joinpath(dirname(@__FILE__), "..", "examples", "differentiation.jl"))
# include(joinpath(dirname(@__FILE__), "..", "examples", "traveling_wave.jl"))
# include(joinpath(dirname(@__FILE__), "..", "examples", "vlasov_evolve.jl"))