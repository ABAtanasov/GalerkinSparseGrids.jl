# We demonstrate the full power of the sparse grid basis
# by successfully evolving a 4+1-dimensional PDE corresponding
# to the phase evolution of some matter distribution. 
# This is the Vlassov-Poisson or "Collisionless Boltzman" equation

# In astrophysical scenarios, this corresponds to 

# The spatial dimension is 2 -> 4D phase space 
# We will use mode order 5 at 5th level resolution along each axis
# Higher resolution at the moment requires a supercomputer 
D = 2; k = 5; n = 5;

# We will evolve from an initial t = 0 to a final t1 = 0.54
t0 = 0; t1 = 0.54

# Initial distribution function in phase space:
fx = x -> sqrt(2*pi) * exp(-2 * pi^2 * norm(x.-1/2)^2)
fv = v -> sqrt(2*pi) * exp(-2 * pi^2 * norm(v)^2)
f0 = x -> fx(2 * (x[1:D] .- 1/2) ) * fv(2 * (x[(D+1):end] .- 1/2) )
fx_modal = vcoeffs_DG(D, k, n, fx)
fv_modal = vcoeffs_DG(D, k, n, fx)
# We use the tensor_construct method as before because of its
# increased accuracy
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

# Finally, we do the full evolution of the Vlasov equation:
vlasov_evolve(D, k, n, 
				m2n, n2p, p2n, n2m, 
				f0_modal, F_point, 
				t0, t1; 
				order="45", points=:specified) 
# points=:specified only saves the first and last set of coefficients 
# in the evolution. This is to save memory.
# If one desires the full set of coefficients, simply remove this keyword arg
