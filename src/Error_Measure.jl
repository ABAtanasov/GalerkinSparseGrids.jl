#------------------------------------------------------
#
# Monte Carlo Integration Techniques for Calculating Error
#
#------------------------------------------------------

#
# More effective than sparse techniques for higher dimensions
#


function monte_carlo(f::Function, D::Int; count = 1000)
    total = 0.0
    for iter in 1:count
		x0 = rand()
		y0 = rand()
        total += f((x0,y0))
	end
    total /= count
end

function monte_carlo2(f::Function, D::Int; batch = 50, Z = 1.0)
    
    means = Float64[]

    while true
        points = Float64[]
        for i in 1:batch
            x0 = rand()
            y0 = rand()
            push!(points, f((x0,y0)))
        end
        if length(means) >1  && (abs(mean(points)-mean(means))/std(means) < Z)
            break
        end
        push!(means, mean(points[end-batch+1:end]))
    end
    return mean(means)
end

function mcerr(f::Function, g::Function, D::Int; count = 1000)
    return sqrt(monte_carlo(x->(f(x)-g(x))^2, D; count=count))
end

function mcerr(f::Function, g::Function, D::Int; batch = 50, Z = 1.0)
	return sqrt(monte_carlo2(x->(f(x)-g(x))^2, D; batch = batch, Z = Z))
end