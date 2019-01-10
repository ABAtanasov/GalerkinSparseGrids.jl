# -----------------------------------------------------------
#
# An example of interpolating a given function in 1D and 2D 
#
# -----------------------------------------------------------

using GalerkinSparseGrids

D = 1
# Random function to interpolate:
truesoln = x -> 1.152 * norm(x)^2 * sin(dot([1.214], x))^3
k_max = 5
n_max = 5

println("Interpolation in ", D, "D.")
println("Going to k_max = ", k_max, ", n_max = ", n_max, ":")

for k_used in 1:k_max
	for n_used in 1:n_max
		# By default, sparse grids are used. Change scheme to "full"  
		# to compare the different accuraces and costs of the 
		# interpolations
		dict = coeffs_DG(D, k_used, n_used, truesoln)
		err = mcerr(x->reconstruct_DG(dict, [x...]), truesoln, D)
		println("(k = ", k_used, ", n = ", n_used, ") : err = ", err)
	end
end

D = 2
# Random function to interpolate:
truesoln = x -> 1.152 * norm(x)^2 * sin(dot([1.214, -0.5], x))^3
k_max = 5
n_max = 5

println("Interpolation in ", D, "D.")
println("Going to k_max = ", k_max, ", n_max = ", n_max, ":")

for k_used in 1:k_max
	for n_used in 1:n_max
		# By default, sparse grids are used. Change scheme to "full"  
		# to compare the different accuraces and costs of the 
		# interpolations
		dict = coeffs_DG(D, k_used, n_used, truesoln)
		err = mcerr(x->reconstruct_DG(dict, [x...]), truesoln, D)
		println("(k = ", k_used, ", n = ", n_used, ") : err = ", err)
	end
end
