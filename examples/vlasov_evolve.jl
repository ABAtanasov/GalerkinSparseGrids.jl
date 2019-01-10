# Constructs the initial 2D Phase space distribution
# Would be ideal to obtain these coefficients using TensorConstruct
# But for D = 2 (4D phase space) this should be fine
function construct_phi(D::Int, k::Int, n::Int, scheme="sparse")
	f = r -> 5/3 * (log(1 + r^2) + 2 * atan(r)/r)
	
	xcoeffs = vcoeffs_DG(D, k, n, x->f(norm(x[1:D])); scheme=scheme)
	pcoeffs = get_one_modal(D, k, n)
	tensor_construct(2*D, k, n, [xcoeffs, pcoeffs])
end



D = 2; k = 5; n = 5;
t0 = 0; t1 = 0.54

# Initial distribution function in phase space:
fx = x -> sqrt(2*pi) * exp(-2 * pi^2 * norm(x.-1/2)^2)
fv = v -> sqrt(2*pi) * exp(-2 * pi^2 * norm(v)^2)
f0 = x -> fx(2 * (x[1:D] .- 1/2) ) * fv(2 * (x[(D+1):end] .- 1/2) )

fx_modal = vcoeffs_DG(D, k, n, fx)
fv_modal = vcoeffs_DG(D, k, n, fx)
f0_modal = tensor_construct(D, k, n, [fx_modal], [fv_modal])

# Constructing the transformation matrices takes the longest time 
# and most memory by far. 
# Once the matrices have been constructed, however, the same 
# matrices can forever be used at a given (D, k, n)
# It would therefore be advantageous to load these in from some
# external file rather than computing them (likely multiple times over)
m2n, n2p = make_modal2point_matrices(D, k, n)
p2n, n2m = make_point2modal_matrices(D, k, n)

# the force as a fucntion of r^2:
# (We use r^2) here because |r| is dicontinuous at the origin and gives
# instabilities for interpolation
F_radial = x -> x == 0 ? 10/9 : 10/3 * (sqrt(x) - atan(sqrt(x)))/(sqrt(x)^3)
fr2 = broadcast(F_radial, get_r2_point(D, k, n, m2n, n2p))
# I have defined F_radial so that multiplying it by x_i 
# gives the ith component of the force vector 
F_point = [get_xi_point(D, i, k, n, m2p) .* fr2 for i in 1:D]

# Finally, we do the full Vlasov evolution:
vlasov_evolve(D, k, n, m2n, n2p, p2n, n2m, f0_modal, F_point, t0, t1; order="45")
