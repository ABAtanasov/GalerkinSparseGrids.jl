# ------------------------------------------------------------------
# During time evolution, we may be concerned that though we do have
# convergence of the function in the L^2 norm, we may have a failure
# of the function's derivative in the sparse basis to correctly 
# converge
# 
# Consequently, we develop a simple filtering method which we can 
# apply at every 
# ------------------------------------------------------------------

# Efficiency criticality: LOW

# Accuracy criticality: MEDIUM

function filter(D::Int, k::Int, n::Int; scheme="sparse", p=5, alpha=15)
    
    cutoff    = get_cutoff(scheme, D, n)
    modes    = ntuple(i-> k, D)
    ls        = ntuple(i->(n+1),D)
    
    V = Float64[]
    for level in CartesianIndices(ls) #This really goes from 0 to l_i for each i
        cutoff(level) && continue
        
        cells = ntuple(i -> 1<<max(0, level[i]-2), D)
        depth = sum
        depth = (1<<sum((level - 1).I))/(1<<(n+1))
        for cell in CartesianIndices(cells)
            for mode in CartesianIndices(modes)
                val = filter(depth; p=p, alpha=alpha)
                push!(V, val)
            end
        end
    end
    return spdiagm(0 => V)
end