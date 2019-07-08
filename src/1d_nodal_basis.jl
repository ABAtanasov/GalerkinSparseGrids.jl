# -----------------------------------------------------------
#
# Implementing a nodal basis for DG sparse grids
#
# -----------------------------------------------------------

#=

We begin by defining analogues of our DG functions
to satisfy the following conditions:

    1) There are a set of collocation points respecting the hierarchical structure
    2) The nodal basis functions at a given level are those that vanish
        on all nodal points at the above levels and vanish on all but one
        of the cells on the level below
    3) Third constraint?

To do this, we _restrict_ to the continuous part of the DG basis, and focus
on piecewise polynomials that vanish on their boundaries. We then pick the
collocation points to be exactly at these boundaries.

For this, it is useful to define lagrange polynomials, and moreover
use k=2, 3, 5.
The k=2 case has our basis is exactly the hat function basis constructed
previously.

For k=2, the coarsest level has nodal points at {0, 1},
    the next level (the first hat) has nodal point at 0.5
    the next level has two cells with points at 0.25, 0.75 etc.

For k=3 the coarsest level has nodal points at {0, 0.5, 1},
    the next level has two nodal functions with points at 0.25, 0.75,
    the next level has two cells and two nodes per cell at all odd multiples of 0.125 etc.

For k=5 the coarsest level has nodal points at {0, 0.25, 0.5, 0.75, 1}
    the next level has four nodal functions with points at all multiples of 0.125
    the next level has two cells and four nodes per cell at all odd multiples of 0.0625 etc.

Notice that k=2,3,5 work because at each cell of the multiresolution basis the space of
    picewise continuous polynomials is k-1 = 2^i for i=0,1,2 in this case.
    A construction for other k would be harder, if not impossible.

=#

# Lagrange polynomials interpolating at the coarsest level
function lag_nodal(k::Int, mode::Int, x::Real)
    (k < 2 || k > 5) && throw(ArgumentError("k is not between 2 and 5"))
    (mode > k || mode < 1) && throw(ArgumentError("mode is not in [1, k]"))

    (x < 0 || x > 1) && return zero(x)
    if k == 2
        mode == 1 && return x
        mode == 2 && return 1 - x
    elseif k == 3
        mode == 1 && return 2 * (x - 1//2) * (x - 1)
        mode == 2 && return -4 * x * (x - 1)
        mode == 3 && return 2 * x * (x - 1//2)
    elseif k == 4
        throw(ArgumentError("k = 4 not implemented"))
    elseif k == 5
        mode == 1 && return 32//3 * (x - 1//4) * (x - 1//2) * (x - 3//4) * (x - 1)
        mode == 2 && return -128//3 * x * (x - 1//2) * (x - 3//4) * (x - 1)
        mode == 3 && return 64 * x * (x - 1//4) * (x - 3//4) * (x - 1)
        mode == 4 && return -128//3 * x * (x - 1//4) * (x - 1//2) * (x - 1)
        mode == 5 && return 32//3 * x * (x - 1//4) * (x - 1//2) * (x - 3//4)
    end
end

# Multiresolution functions for the levels below
function h_nodal(k::Int, mode::Int, x::Real)
    (k < 2 || k > 5) && throw(ArgumentError("k is not between 2 and 5"))
    (mode > k || mode < 1) && throw(ArgumentError("mode is not in [1, k]"))

    if k == 2
        mode == 1 && return max(zero(x), 1-abs(x))
        if mode == 2
            (abs(x) > 1 || abs(x) == 0) && return zero(x)
            x < 0 && return (x + 1)
            x > 0 && return zero(x)# (x - 1)
        end
    elseif k == 3
        mode == 1 && return max(zero(x), -4 * x * (x + 1))
        mode == 2 && return max(zero(x), -4 * x * (x - 1))
        if mode == 3
            (abs(x) > 1 || abs(x) == 0) && return zero(x)
            x < 0 && return 2 * (x + 1) * (x + 1//2)
            x > 0 && return zero(x) # -2 * (x - 1) * (x - 1//2)
        end
    elseif k == 4
        throw(MethodError("k = 4 not implemented"))
    elseif k == 5
        mode == 1 && return (x < -1 || x > 0) ? zero(x) : -128//3 * x * (x + 1) * (x + 1//2) * (x + 1//4)
        mode == 2 && return (x < -1 || x > 0) ? zero(x) : -128//3 * x * (x + 1) * (x + 1//2) * (x + 3//4)
        mode == 3 && return (x > 1  || x < 0) ? zero(x) : -128//3 * x * (x - 1) * (x - 1//2) * (x - 3//4)
        mode == 4 && return (x > 1  || x < 0) ? zero(x) : -128//3 * x * (x - 1) * (x - 1//2) * (x - 1//4)
        if mode == 5
            (abs(x) > 1 || abs(x) == 0) && return zero(x)
            x < 0 && return 32//3 * (x + 1) * (x + 3//4) * (x + 1//2) * (x + 1//4)
            x > 0 && return zero(x) #32//3 * (x - 1) * (x - 3//4) * (x - 1//2) * (x - 1//4)
        end
    end
end

# This yields our nodal basis:
function v_nodal(k::Int, level::Int, cell::Int, mode::Int, x::Real)
    if level==0
        return lag_nodal(k, mode, x)
    else
        return h_nodal(k, mode, (1<<level)*x - (2*cell-1))
    end
end

function v_nodal(k::Int, level::Int, cell::Int, mode::Int)
    return x->v_nodal(k, level, cell, mode, x)
end

# -----------------------------------------------------------
# Made for evaluating the dicontinuous elements in the nodal
# basis at their points of discontinuity
# -----------------------------------------------------------

# Evaluates the modal basis at a given point, averaging over jumps
function eval_v_periodic(k::Int, level::Int, cell::Int, mode::Int, x::T) where {T<:Real}
    xl = max(x - 1e-16, 0)
    xr = min(x + 1e-16, 1)
    (xl == 0) && (xl = 1 - 1e-16)
    (xr == 1) && (xr = 1e-16)
    return 1//2 * (v(k, level, cell, mode, xl) + v(k, level, cell, mode, xr))
end

# Find the size of the discontinuity of the modal basis at a given point
function gap_v(k::Int, level::Int, cell::Int, mode::Int, x::T) where {T<:Real}
    xl = max(x - 1e-16, 0)
    xr = min(x + 1e-16, 1)
    return 1//2 * (v(k, level, cell, mode, xr) - v(k, level, cell, mode, xl))
end

# Evaluates the modal basis at a given point, averaging over jumps
function eval_v_nodal_periodic(k::Int, level::Int, cell::Int, mode::Int, x::T) where {T<:Real}
    xl = max(x - 1e-30, 0)
    xr = min(x + 1e-30, 1)
    (xl == 0) && (xl = 1 - 1e-30)
    (xr == 1) && (xr = 1e-30)
    return 1//2 * (v_nodal(k, level, cell, mode, xl) + v_nodal(k, level, cell, mode, xr))
end

# Find the size of the discontinuity of the modal basis at a given point
function avg_v_nodal(k::Int, level::Int, cell::Int, mode::Int, x::T) where {T<:Real}
    xl = x - 1e-16
    xr = x + 1e-16
    (xl < 0) && (xl = xr)
    (xr > 1) && (xr = xl)
    return 1//2 * (v_nodal(k, level, cell, mode, xr) + v_nodal(k, level, cell, mode, xl))
end

# Find the size of the discontinuity of the modal basis at a given point
function gap_v_nodal(k::Int, level::Int, cell::Int, mode::Int, x::T) where {T<:Real}
    xl = x - 1e-16
    xr = x + 1e-16
    (xl < 0) && (xl = xr)
    (xr > 1) && (xr = xl)
    return 1//2 * (v_nodal(k, level, cell, mode, xr) - v_nodal(k, level, cell, mode, xl))
end

function eval_v_nodal_right(k::Int, level::Int, cell::Int, mode::Int, x::T) where {T<:Real}
    xr = x + 1e-16
    return v_nodal(k, level, cell, mode, xr)
end

function eval_v_nodal_left(k::Int, level::Int, cell::Int, mode::Int, x::T) where {T<:Real}
    xl = x - 1e-16
    return v_nodal(k, level, cell, mode, xl)
end

# -----------------------------------------------------------
# Next we define the sparse matrices transforming between
# the heirarchical, nodal, and point bases
# -----------------------------------------------------------

# Huge TODO:
#     Replace all of these for loops and multiple function calls
#     by cleverly using the iterable interface in julia

# Evaluation of a given nodal function along all relevant gridpoints
function eval_points_1D(k::Int, max_level::Int, level::Int, cell::Int, mode::Int)
    i = 1
    I = Int[]
    V = Float64[]
    for l in 0:max_level
        l == 0 ? (node_range = 0:(1//(k-1)):(1 - 1//(k-1))) :
                 (node_range = 1//(2*(k-1)):1//(k-1):(1-1//(2*(k-1))))

        for c in 1:1<<max(0, l-1)
            for node in node_range
                x = (c - 1 + node)//(1<<max(0, l-1))
                val = eval_v_nodal_right(k, level, cell, mode, x)
                (val != zero(x)) && (push!(I, i); push!(V, val))
                i += 1
            end

            if l == 0
                x = 1
                val = eval_v_nodal_left(k, level, cell, mode, x)
            else
                x = (c - 1//2)/(1<<max(0, l-1))
                val = eval_v_nodal_left(k, level, cell, mode, x)
            end
            (val != zero(x)) && (push!(I, i); push!(V, val))
            i += 1
        end
    end
    return sparsevec(I, V)
end

# Sparse matrix converting nodal basis element to its value
# on the collocation points
function nodal2points_1D(k::Int, max_level::Int; atol=1e-15)
    I = Int[]
    J = Int[]
    V = Float64[]
    j = 1

    for level in 0:max_level
        mode_range = 1:k

        for cell in 1:1<<max(0, level-1)
            for mode in mode_range
                coeffs = eval_points_1D(k, max_level, level, cell, mode)
                for (i, val) in zip(findnz(coeffs)...)
                    push!(I, i)
                    push!(J, j)
                    push!(V, val)
                end
                j += 1
            end
        end
    end
    return threshold(sparse(I, J, V), atol)
end

# Sparse matrix inverting the above construction
# Converting a set of values on the collocation points to the DG space
# in terms of the nodal basis
function points2nodal_1D(k::Int, max_level::Int; atol=1e-15)
    n2p = nodal2points_1D(k, max_level)
    return threshold(inv(Matrix(n2p)), atol)
end

# Converts a given nodal basis element to a sum in the standard
# hieararchical DG scheme
# TODO: Remove hquadrature and do this directly
function nodal2modal_1D(k::Int, level::Int, cell::Int, mode::Int; rtol=1e-10, atol=1e-15, max_evals=1500)
    I = Int[]
    V = Float64[]
    i = 1

    for hier_level in 0:level
        hier_cells = 1<<max(0, hier_level-1)
        for hier_cell in 1:hier_cells
            for hier_mode in 1:k
                x_min = (hier_cell-1)//hier_cells
                x_max = (hier_cell)//hier_cells
                val = hquadrature(x->v_nodal(k, level, cell, mode, x)*
                                    v(k, hier_level, hier_cell, hier_mode, x),
                                    x_min, x_max,
                                    reltol=rtol, abstol=atol, maxevals=max_evals)[1]
                (abs(val) > eps(atol)) && (push!(I, i); push!(V, val))
                i += 1
            end
        end
    end
    return dropzeros!(sparsevec(I, V))
end

function nodal2modal_1D(k::Int, max_level::Int; rtol=1e-10, atol=1e-15, max_evals=1500)
    I = Int[]
    J = Int[]
    V = Float64[]
    j = 1

    for nodal_level in 0:max_level
        num_nodes = k
        nodal_cells = 1<<max(0, nodal_level-1)
        for nodal_cell in 1:nodal_cells
            for nodal_mode in 1:num_nodes
                coeffs = nodal2modal_1D(k, nodal_level, nodal_cell, nodal_mode;
                                        rtol=rtol,
                                        atol=atol,
                                        max_evals=max_evals)
                for (i, val) in zip(findnz(coeffs)...)
                    push!(I, i)
                    push!(J, j)
                    push!(V, val)
                end
                j += 1
            end
        end
    end
    return threshold(sparse(I, J, V), atol)
end

# Accounting for ambiguity in the definition of the hierarchical
# basis across regions
function eval_v(k, level, cell, mode, x)
    xl = max(x - 1e-16, 0)
    xr = min(x + 1e-16, 1)
    (xl == 0) && (xl = 1 - 1e-16)
    (xr == 1) && (xr = 1e-16)
    return 0.5 * (v(k, level, cell, mode, xl) + v(k, level, cell, mode, xr))
end

function nodal2pos_1D(k, max_level, nodal_level, nodal_cell, nodal_mode;
                   rtol=1e-10, atol=1e-15, max_evals=1500)
    I = Int[]
    V = Float64[]
    i = 1

    nodal_min = (nodal_cell-1)/(1<<max(0, nodal_level-1))
    nodal_med = (nodal_cell-1//2)/(1<<max(0, nodal_level-1))
    nodal_max = (nodal_cell)/(1<<max(0, nodal_level-1))
    for c in 1:1<<max_level
        pos_min = (c-1)/(1<<max_level)
        pos_med = (c-1//2)/(1<<max_level)
        pos_max = c/(1<<max_level)

        sm = 1e-20

        if (pos_min > nodal_max) || (pos_max < nodal_min)
            i += k
            continue
        end
        for m in 1:k
            val1  = hquadrature(x->v_nodal(k, nodal_level, nodal_cell, nodal_mode, x)*basis(max_level, c, m, x),
                              pos_min + sm, pos_med - sm,
                              rtol=rtol, atol=atol, maxevals=max_evals)[1]

            val2 = hquadrature(x->v_nodal(k, nodal_level, nodal_cell, nodal_mode, x)*basis(max_level, c, m, x),
                              pos_med + sm, pos_max - sm,
                              rtol=rtol, atol=atol, maxevals=max_evals)[1]

            val = val1 + val2
            (abs(val) > eps(atol)) && (push!(I, i); push!(V, val))
            i += 1
        end
    end
    return sparsevec(I, V)
end

function nodal2pos_1D(k, max_level; rtol=1e-10, atol=1e-15, max_evals=1500)
    I = Int[]
    J = Int[]
    V = Float64[]
    j = 1

    for nodal_level in 0:max_level
        num_nodes = k
        nodal_cells = 1<<max(0, nodal_level-1)
        for nodal_cell in 1:nodal_cells
            for nodal_mode in 1:num_nodes
                coeffs = nodal2pos_1D(k, max_level, nodal_level, nodal_cell, nodal_mode;
                                        rtol=rtol,
                                        atol=atol,
                                        max_evals=max_evals)
                for (i, val) in zip(findnz(coeffs)...)
                    push!(I, i)
                    push!(J, j)
                    push!(V, val)
                end
                j += 1
            end
        end
    end
    return threshold(sparse(I, J, V), atol)
end

function pos2nodal_1D(k, max_level; rtol=1e-10, atol=1e-12, max_evals=1500)
    n2pos = nodal2pos_1D(k, max_level; rtol=rtol, atol=atol, max_evals=max_evals)
    return threshold(inv(Matrix(n2pos)), atol)
end


# Final function, combining all of the above:
function transform_1D(k::Int, n::Int, from::String, to::String; atol=1e-12)
    if from == "nodal" && to == "pos"
        return nodal2pos_1D(k, n; atol=atol)
    elseif from == "pos" && to == "nodal"
        return pos2nodal_1D(k, n; atol=atol)
    elseif from == "nodal" && to == "points"
        return nodal2points_1D(k, n; atol=atol)
    elseif from == "points" && to == "nodal"
        return points2nodal_1D(k, n; atol=atol)
    elseif from == "nodal" && to == "modal"
        return threshold(pos2hier(k, n) * nodal2pos_1D(k, n), atol)
    elseif from == "modal" && to == "nodal"
        return threshold(pos2nodal_1D(k, n) * hier2pos(k, n), atol)
    elseif from == "modal" && to == "pos"
        return hier2pos(k, n; atol=atol)
    elseif from == "pos" && to == "modal"
        return pos2hier(k, n; atol=atol)
    elseif from == "points" && to == "pos"
        return threshold(nodal2pos_1D(k, n) * points2nodal_1D(k, n), atol)
    elseif from == "pos" && to == "points"
        return threshold(nodal2points_1D(k, n) * pos2nodal_1D(k, n), atol)
    elseif from == "modal" && to == "points"
        m2n = transform_1D(k, n, "modal", "nodal"; atol=atol)
        n2p = transform_1D(k, n, "nodal", "points"; atol=atol)
        return threshold(n2p * m2n, atol)
    elseif from == "points" && to == "modal"
        p2n = transform_1D(k, n, "points", "nodal"; atol=atol)
        n2m = transform_1D(k, n, "nodal", "modal"; atol=atol)
        return threshold(n2m * p2n, atol)
    else throw(ArgumentError("A basis was undefined"))
    end
end
