# -----------------------------------------------------
# This script gives the methods for maniupulating the
# Cartesian indices corresponding to multilevels that
# are either kept or cut in various schemes of
# "sparsification"
# -----------------------------------------------------

# Efficiency Criticality: LOW
# This is not evaluated very often

# Accuracy Criticality: N/A
# No float manipulation



# Gives the appropriate boolean cutoff corresponding
# to a given scheme, e.g. sparse basis, full basis,
# and in the future possibly the energy basis of
# Bungartz and Griebel

function cutoff(::Val{:sparse}, x::CartesianIndex{D}, n::Int) where {D}
    sum(x.I) > n+D
end

function cutoff(::Val{:full}, x::CartesianIndex{D}, n::Int) where {D}
    false
end
