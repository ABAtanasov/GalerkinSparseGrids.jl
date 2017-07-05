using GalerkinSparseGrids
using Base.Test
using Cubature
using ODE

tests = [
	"elementary.jl",
	"hat_basis.jl",
	"hier_DG.jl",
	"vhier_DG.jl",
	"differentiation.jl",
	"solvers.jl"
]

for filename in tests
	name = first(splitext(filename))
	@testset "$name" begin
		include(filename)
	end
end