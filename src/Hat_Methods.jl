#------------------------------------------------------
#
# All methods involving the naive non-galerkin hat basis
# of Griebel et al. 
#
#------------------------------------------------------

# Efficiency criticality: MEDIUM
# Important for calculating error for testing package
# Not critical for user
# Potential use to replace multidimensional integration 
# from hcubature with more efficient method


#------------------------------------------------------
# Defining the Hat Functions
#------------------------------------------------------

#Our original function which we will shift and scale
function hat(x::Real)
    b=abs(x)
    if b>1
        return zero(x)
    end
    return 1-b
end

#Our 1-D basis
function hat(l::Int,i::Int,x::Real)
    if l==-1				#The constant function (doesn't vanish on boundary)
        return one(x)		
    end
    if l==0					#The linear function (vanishes on only one boundary)
        return x
    end
    return hat((1<<l)*x - (2i - 1)) #Shifts and scalings of the original hat function
end

#Tensor product construction: our d-D basis
function hat{D,T<:Real}(ls::NTuple{D,Int},is::NTuple{D,Int},xs::NTuple{D,T}) 
    ans=one(eltype(xs))
    for k = 1:length(ls)
        ans *= hat(ls[k],is[k],xs[k])
    end
    return ans
end

#Lagrange function that spikes up only at a point
function delta(l::Int,i::Int,x::Float64)
    return hat((1<<l) * x - i)
end

# ψ{D}(ls::NTuple{D,Int}, is::NTuple{D,Int}, ....)
#Tensor product construction of Lagrange basis
function delta{D,T<:Real}(ls::NTuple{D,Int}, is::NTuple{D,Int}, xs::NTuple{D,T})
    ans=one(eltype(xs))
    for k = 1:length(ls)
        ans *= delta(ls[k],is[k],xs[k])
    end
    return ans
end


#------------------------------------------------------
# Lagrange Basis
#------------------------------------------------------
function standard_coefficients{D}(f::Function, ls::NTuple{D,Int})
    positions = ntuple(i -> (1<<ls[i])+1,D)
    coeffs = zeros(Float64, positions)
    for place in CartesianRange(positions)
        x = ntuple(i-> (2.0^-ls[i])*(place[i]-1),D)
        coeffs[place] = f(x) 
    end
    return coeffs
end

function standard_reconstruct{D,T<:Real}(coefficients::AbstractArray, ls::NTuple{D,Int}, xs::NTuple{D,T})
    positions = ntuple(i -> (1<<ls[i])+1,D)
    value=0.0
    for place in CartesianRange(positions)
        is = ntuple(i->place[i]-1,D)
        value += coefficients[place]*delta(ls, is, xs)
    end
    return value
end


#------------------------------------------------------
# Methods for obtaining position (lagrange basis) 
# coefficients
#------------------------------------------------------

# Given a level and place describing a sub-interval 
# gives the x-coordinate for the center of that sub-interval. 
# This is used for obtaining the coefficient corresponding to 
# a specific level and place in a position (lagrange basis).
# This function's implementation is multi-dimensional
function get_position{D}(level::CartesianIndex{D},place::CartesianIndex{D})
    xs = [(0.5)^pos(level[i]-2) *(2*place[i]-1) for i in 1:D]
    return xs
end

# Takes an array xs of length n and makes a new n+1 length array
# obtained from joining j to the beginning of array xs
# without altering xs
function form_array{T<:Real}(j::Real,xs::Array{T})
    return [i==1?j:xs[i-1] for i in 1:(length(xs)+1)]
end

# Takes a function f of an n-D vector, and returns a function f' of an (n-1)-D
# vector x', so that f'(x') = f( [j, x']). That is, the first coordinate of x has been
# fixed to equal j
function project_function(f::Function, coordinate::Real)
    return (x -> f(form_array(coordinate,x)))
end

# This is how we get the coefficients. 
# This recursive implementation is the easiest way to write it. 
# 
# TODO work on efficiency improvements here. This will be the basis for sparse integration
#
function get_coefficient{D,T<:Real}(f::Function, 
									level::NTuple{D,Int}, 
									x::Array{T}) 
									
    if D==1					#This is the base case for 1D on the interval [0,1]
        if level[1]==1
            return f([0.])
        elseif level[1]==2
            return f([1.])-f([0.])
        else
            middle = f(x)
            x[1] += (0.5)^(level[1]-2)
            after = f(x)
            x[1] -= 2*(0.5)^(level[1]-2)
            before = f(x)
            x[1] += (0.5)^(level[1]-2) #This makes sure x is unchanged in the end
            return middle - 0.5 * before - 0.5 * after
        end
    else					# Otherwise, we will take appropriate differences of differences
							# depending on the first coefficient (treating that as our first dim)
        if level[1]==1
            return get_coefficient(project_function(f,0), level[2:end], x[2:end])
        elseif level[1]==2
            after  = get_coefficient(project_function(f,1), level[2:end], x[2:end])
            before = get_coefficient(project_function(f,0), level[2:end], x[2:end])
            return after - before
        else
            middle = get_coefficient(project_function(f,x[1]), level[2:end], x[2:end])
            x[1] += (0.5)^(level[1]-2)
            after  = get_coefficient(project_function(f,x[1]), level[2:end], x[2:end])
            x[1] -= 2*(0.5)^(level[1]-2)
            before = get_coefficient(project_function(f,x[1]), level[2:end], x[2:end])
            x[1] += (0.5)^(level[1]-2) #This makes sure x is unchanged in the end
            return middle - 0.5 * before - 0.5 * after
        end
    end
end

#------------------------------------------------------
# Hierarchical Coefficients in n-D
#------------------------------------------------------
function hier_coefficients{D}(f::Function, ls::NTuple{D,Int})
    coeffs = Dict{CartesianIndex{D}, Array{Float64,D}}()
	# We will make a dictionary, that given a level (represented by a cartesian index), 
	# will lead to a list of coefficients representing the places at that level
    for level in CartesianRange(ls)     			# We go from from -1 to l_i for each i,
													# This gives a shift by 2
        ks = ntuple(i -> 1<<pos(level[i]-3), D) 	# Dimension of the ks
        level_coeffs = zeros(ks)		#the coefficients for all the places at this level
        for place in CartesianRange(ks)
            x = get_position(level, place)  # This is the position array corresponding to a place.
											# It will get modified in calculating the coefficient
            level_coeffs[place]=get_coefficient(f,ntuple(i -> level[i], D),x)
        end
        coeffs[level] = level_coeffs
    end
    return coeffs
end


#------------------------------------------------------
# Sparse Coefficients in n-D
#------------------------------------------------------
function sparse_coefficients(f::Function, n::Int, D::Int)
    coeffs = Dict{CartesianIndex{D}, Array{Float64,D}}()
    ls = ntuple(i -> (n+2), D)			# n + 2 because we go From -1 to l_i for each i
    for level in CartesianRange(ls) 	
        diag_level = 0;
        for i in 1:D
            diag_level += level[i]
        end
        if diag_level > n + 2*D 		# If we're past the relevant levels
            continue					# don't compute coeffs
        else  							#Otherwise we'll go ahead and compute them
            ks = ntuple(i -> 1<<pos(level[i]-3), D)  
            level_coeffs = zeros(Float64, ks)
            for place in CartesianRange(ks)
                x = get_position(level, place) 
                level_coeffs[place]=get_coefficient(f,ntuple(i -> level[i], D),x)
            end
            coeffs[level] = level_coeffs
        end
    end
    return coeffs
end

#------------------------------------------------------
# Reconstructing from Coefficients in n-D
#------------------------------------------------------
#this will work for both sparse and full cases
#it also generalizes decently nicely to the DG case
function reconstruct{D,T<:Real}(coefficients::Dict{CartesianIndex{D},
		 						Array{Float64,D}}, x::NTuple{D,T})
    value = 0.0
    for key in keys(coefficients)	#For every level that has coefficients
        level = ntuple(i->key[i]-2,D)	# Get the actual level corresponding to that CartesianIndex
        place = ntuple(i->hat_index(x[i],level[i]),D)  # Get the relevant place for our position x
        value += (coefficients[key])[CartesianIndex{D}(place)]*hat(level,place,x) 
		#get the appropriate coefficient and evaluate the appropriate ϕ at x
    end
    return value	#return the sum of all the relevant hat functions at that place x
end