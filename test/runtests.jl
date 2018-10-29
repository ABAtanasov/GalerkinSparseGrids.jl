using GalerkinSparseGrids
using Test
using HCubature
using ODE
using SparseArrays

tests = [
	"elementary.jl",
	"hier_DG.jl",
	"vhier_DG.jl",
	"differentiation.jl",
	"solvers.jl"
]
@testset "GalerkinSparseGrids" begin
	for filename in tests
		name = first(splitext(filename))
			include(filename)
	end
end
