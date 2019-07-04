# We demonstrate the full power of the sparse grid basis
# by successfully evolving a 4+1-dimensional PDE corresponding
# to the phase evolution of some matter distribution.
# This is the Vlassov-Poisson or "Collisionless Boltzman" equation
using GalerkinSparseGrids
using LinearAlgebra
using ODE

# The spatial dimension is 2, meaning a 4D phase space
# We will use mode order 5 at 5th level resolution along each axis
# Higher resolution at the moment requires a supercomputer
D = 2; k = 5; n = 5
# D = 2; k = 5; n = 2

# We will evolve from an initial t = 0 to a final t1 = 0.54
t0 = 0.0; t1 = 0.54; nout = 101
# t0 = 0.0; t1 = 0.0054; nout = 2

# Initial distribution function in phase space:
@inline fx(x) = exp(-2 * pi^2 * norm(x.-1/2)^2)
@inline fv(v) = exp(-2 * pi^2 * norm(v)^2)
fx_modal = coeffs_DG(1, k, n, fx)
fv_modal = coeffs_DG(1, k, n, fv)
# We use the tensor_construct method as before because of its
# increased accuracy
f0_modal_array = [i <= D ? fx_modal : fv_modal for i in 1:2*D]
f0_modal = 2*pi*D2V(2*D,k,n,tensor_construct(2*D, k, n, f0_modal_array))

# Constructing the transformation matrices takes the longest time
# and most memory by far.
# Once the matrices have been constructed, however, the same
# matrices can forever be used at a given (D, k, n)
# It would therefore be advantageous to load these in from some
# external file rather than computing them (likely multiple times over)
m2n, n2p = make_modal2point_matrices(2*D, k, n)
p2n, n2m = make_point2modal_matrices(2*D, k, n)

# the force as a function of r^2:
# (We use r^2) here because |r| is dicontinuous at the origin and gives
# instabilities for interpolation
F_radial = x -> x == 0 ? 10/9 : 10/3 * (sqrt(x) - atan(sqrt(x)))/(sqrt(x)^3)
fr2 = broadcast(F_radial, get_r2_point(2*D, k, n, m2n, n2p))
# I have defined F_radial so that multiplying it by x_i
# gives the ith component of the force vector
F_point = [get_xi_point(2*D, i, k, n, m2n, n2p) .* fr2 for i in 1:D]

# Finally, we do the full evolution of the Vlasov equation:
vlasov_evolve(D, k, n,
                m2n, n2p, p2n, n2m,
                f0_modal, F_point,
                t0, t1, nout;
                order="45", points=:specified)
# points=:specified only saves the first and last set of coefficients
# in the evolution. This is to save memory.
# If one desires the full set of coefficients, simply remove this keyword arg
