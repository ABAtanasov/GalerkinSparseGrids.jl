# -----------------------------------------------------------
#
# An example of contructing the initial data for a traveling
# wave in D-dimensions using the tensor_construct method
# and the wave_evolve in PDEs.jl
#
# -----------------------------------------------------------

using GalerkinSparseGrids
using LinearAlgebra

# -----------------------------------------------------
# Generates the coefficients for a waveform
# A cos(\vec{k} \cdot \vec{x} + \phi), \vec{k} = 2 \pi \vec{m}
# using an interpolation of type (k,n)
# and periodic boundary conditions in D-dimensionss
# -----------------------------------------------------
function cos_coeffs(k::Int, n::Int, m::AbstractArray{T,1};
                    scheme="sparse", phase = 0.0, A = 1.0) where T
    D = length(m)
    wavenumber = 2*pi*m

    # Begin with writing cos(\sum_i k_i x_i) as a sum of products
    # of sines/cosines of individual k_j x_j
    sines   = [i==1 ? x->sin(wavenumber[i]*x[1]+phase) : x->sin(wavenumber[i]*x[1]) for i in 1:D]
    cosines = [i==1 ? x->cos(wavenumber[i]*x[1]+phase) : x->cos(wavenumber[i]*x[1]) for i in 1:D]

    sine_dicts   = [coeffs_DG(1, k, n, sines[i]) for i in 1:D]
    cosine_dicts = [coeffs_DG(1, k, n, cosines[i]) for i in 1:D]

    ansVect = zeros(get_size(Val(D), k, n, Val(Symbol(scheme))))

    for SCs in CartesianIndices(ntuple(q->2, D))
        num_sines = sum([SCs[i]-1 for i in 1:D])
        if num_sines % 2 == 1
            continue
        end
        sign = num_sines%4==0 ? 1 : -1

        coeff_array = [SCs[i]==1 ? cosine_dicts[i] : sine_dicts[i] for i in 1:D]
        productDict = tensor_construct(D, k, n, coeff_array; scheme=scheme)
        productVect = D2V(D, k, n, productDict; scheme=scheme)

        ansVect += sign * productVect
    end

    return A * ansVect
end

# -----------------------------------------------------
# The same as above, but using sin
# -----------------------------------------------------
function sin_coeffs(k::Int, n::Int, m::Array{Int,1};
                    scheme="sparse", phase=0.0, A=1.0)
    return cos_coeffs(k, n, m; scheme=scheme, phase=phase-pi/2, A=A)
end

# -----------------------------------------------------
# Returns the data for a traveling wave
# using the above methods
# ----------------------------------------------------
function traveling_wave(k::Int, n::Int, m::Array{Int,1};
                        scheme="sparse", phase=0.0, A=1.0)
    wavenumber = 2*pi*m
    frequency = sqrt(dot(wavenumber,wavenumber))

    # u(x) = A * cos(dot(wavenumber,x) + phase)
    # v(v) = A * frequency * sin(dot(k,x) + phase)
    u0_coeffs = cos_coeffs(k, n, m; scheme=scheme, phase=phase, A=A)
    v0_coeffs = sin_coeffs(k, n, m; scheme=scheme, phase=phase, A=A*frequency)
    return (u0_coeffs, v0_coeffs)
end

# -----------------------------------------------------
# Main routine:
# -----------------------------------------------------

# Modify m to change the wavenumber.
# We use integer entries here because of periodic boundary conditions
m = [1,2,-1]
truesoln = x -> cos(2*pi*(dot(m,x) - sqrt(dot(m,m))*0.54))
k_max = 5
n_max = 6
D = length(m)

t0 = 0;
t1 = 0.54;

println("Wave Evolution in ", D, "D.")
println("Going to k_max = ", k_max, ", n_max = ", n_max, ":")

for k_used in 1:k_max
    for n_used in 1:n_max
        f0coeffs, v0coeffs = traveling_wave(k_used, n_used, m)
        soln = wave_evolve(D, k_used, n_used, f0coeffs, v0coeffs, t0, t1)
        dict = V2D(D, k_used, n_used, soln[2][end])

        err = mcerr(x->reconstruct_DG(dict, [x...]), truesoln, D)
        println("(k = ", k_used, ", n = ", n_used, ") : err = ", err)
    end
end


