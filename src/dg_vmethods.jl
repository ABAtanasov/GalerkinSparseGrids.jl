# ------------------------------------------------------------------
# It is advantageous to have the coeffs as
# one single vector in order to make use of BLAS
# and related libraries when defining operators on
# our space of functions
#
# In this script, we "vectorize" our functions,
# not in the regular meaning of the word, but in the
# sense that dictionaries will all be recelld by vectors
#
# The main difficulty with this is that, unlike with dictionaries
# there is no easy way to go from (level, cell, mode) to a
# corresponding index in a 1-D vector.
#
# For this, we have the reference functions:
# full/sparse_referenceV2D/D2V
# D2V generates a dict that, upon input of a level, cell, mode
# gives the corresponding index in a vector
# V2D generates a vector, with row i having the three numbers
# level, cell, mode corresponding to index i in the vector
#
# These methods work in all dimensions, and there are ones for both
# full and sparse grids.
#
# The entire script culminates in a final result: a matrix
# representation of the derivative operator
# (both in full and sparse bases)
# ------------------------------------------------------------------


# Calculates the exact dimension of interpolating
# basis functions in the full or sparse grid schemes
# using the degree k-1 Galerkin polynomials
# with multi-resolution up to n
function get_size(::Val{D}, k::Int, n::Int, scheme::Val{Scheme}=Val(:sparse)) where {D, Scheme}
    ls      = ntuple(i->(n+1), D)
    size    = 0
    for level in CartesianIndices(ls)
        cutoff(scheme, level, n) && continue

        ks = ntuple(q -> 1<<max(0, level[q]-2), D)
        size += prod(ks)*k^D
    end
    return size
end


function D2V(d::Int, k::Int, n::Int, coeffs::Dict{CartesianIndex{D}, <:AbstractArray{<:AbstractArray{T, D}, }};
    scheme="sparse") where {D, T <: Real}
    @assert d == D
    D2V(Val{D}(), k, n, coeffs, Val(Symbol(scheme)))
end
function D2V(::Val{D}, k::Int, n::Int, coeffs::Dict{CartesianIndex{D}, <:AbstractArray{<:AbstractArray{T, D}, }},
    scheme::Val{Scheme}) where {D, T <: Real, Scheme}

    size        = get_size(Val{D}(), k, n, scheme)
    vect        = Array{T}(undef, size)
    modes       = ntuple(i-> k, D)
    ls            = ntuple(i->(n+1), D)
    j = 1
    for level in CartesianIndices(ls) #This really goes from 0 to l_i for each i
        cutoff(scheme, level, n) && continue

        ks = ntuple(q -> 1<<max(0, level[q]-2), D) #This sets up a specific k+1 vector
        for cell in CartesianIndices(ks)
            for mode in CartesianIndices(modes)
                vect[j] = coeffs[level][cell][mode]
                j += 1
            end
        end
    end
    return vect
end


function V2D(D::Int, k::Int, n::Int, vect::Array{T}; scheme="sparse") where {T <: Real}
    V2D(Val{D}(), k, n, vect, Val(Symbol(scheme)))
end
function V2D(::Val{D}, k::Int, n::Int, vect::Array{T}, scheme::Val{Scheme}) where {D, T <: Real, Scheme}
    coeffs        = Dict{CartesianIndex{D}, Array{Array{T,D},D}}()
    modes        = ntuple(q-> k, D)
    ls            = ntuple(i->(n+1), D)
    j = 1
    for level in CartesianIndices(ls) #This really goes from 0 to l_i for each i
        cutoff(scheme, level, n) && continue

        ks = ntuple(q -> 1<<max(0, level[q]-2), D)  #This sets up a specific k+1 vector
        level_coeffs = Array{Array{T}}(undef, ks) #all the coefficients at this level
        for cell in CartesianIndices(ks)
            cell_coeffs = Array{T}(undef, modes)
            for mode in CartesianIndices(modes)
                cell_coeffs[mode] = vect[j]
                j += 1
            end
            level_coeffs[cell] = cell_coeffs
        end
        coeffs[level] = level_coeffs
    end
    return coeffs
end

function D2Vref(::Val{D}, k::Int, n::Int, scheme::Val{Scheme}) where {D, Scheme}
    size        = get_size(Val{D}(), k, n, scheme)
    dict        = Dict{NTuple{3,CartesianIndex{D}}, Int}()
    modes        = ntuple(q-> k, D)
    ls            = ntuple(i->(n+1), D)
    j = 1
    for level in CartesianIndices(ls)
        cutoff(scheme, level, n) && continue

        ks = ntuple(q -> 1<<max(0, level[q]-2), D)  #This sets up a specific k+1 vector
        lvl = ntuple(i -> level[i]-1,D)
        for cell in CartesianIndices(ks)
            for mode in CartesianIndices(modes)
                dict[(level, cell, mode)] = j
                j += 1
            end
        end
    end
    return dict
end

function V2Dref(::Val{D}, k::Int, n::Int, scheme::Val{Scheme}) where {D, Scheme}
    size        = get_size(Val{D}(), k, n, scheme)
    vect        = Array{NTuple{3,CartesianIndex{D}}}(undef, size)
    modes       = ntuple(q-> k, D)
    ls            = ntuple(i->(n+1), D)
    j = 1
    for level in CartesianIndices(ls)
        cutoff(scheme, level, n) && continue

        ks = ntuple(q -> 1<<max(0, level[q]-2), D)  #This sets up a specific k+1 vector
        lvl = ntuple(i -> level[i]-1,D)
        for cell in CartesianIndices(ks)
            for mode in CartesianIndices(modes)
                vect[j] = (level, cell, mode)
                j += 1
            end
        end
    end
    return vect
end


# -----------------------------------------------------
# Let's now make the coefficient operators work on
# and return vectors
# -----------------------------------------------------
function vcoeffs_DG(D::Int, k::Int, n::Int, f::Function;
                                rtol=REL_TOL, atol=ABS_TOL,
                                maxevals=MAX_EVALS,
                                scheme="sparse")
    vcoeffs_DG(Val(D), k, n, f, rtol, atol, maxevals, Val(Symbol(scheme)))
end
function vcoeffs_DG(::Val{D}, k::Int, n::Int, f::Function,
                                rtol, atol,
                                maxevals,
                                scheme::Val{Scheme}) where {D, Scheme}
    len            = get_size(Val{D}(), k, n, scheme)
    coeffs        = Array{Float64}(undef, len)
    modes        = ntuple(i-> k, D)
    ls            = ntuple(i-> (n+1), D)
    j = 1
    for level in CartesianIndices(ls)  # This really goes from 0 to l_i for each i,
        cutoff(scheme, level, n) && continue

        ks = ntuple(i -> 1<<max(0, level[i]-2), D)  #This sets up a specific k+1 vector
        lvl = ntuple(i -> level[i]-1,D)
        for cell in CartesianIndices(ks)
            for mode in CartesianIndices(modes)
                coeffs[j] = get_coefficient_DG(k, lvl, cell, mode, f;
                                                rtol=rtol, atol=atol,
                                                maxevals=maxevals)
                j += 1
            end
        end
    end
    return coeffs
end
