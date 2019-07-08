# -----------------------------------------------------
#
# Monte Carlo Integration Techniques for Calculating Error
#
# -----------------------------------------------------

#
# More effective than sparse techniques for higher dimensions
#


function monte_carlo(f::Function, D::Int; count = 1000)
    total = 0.0
    for iter in 1:count
        x = ntuple(i -> rand(), D)
        total += f(x)
    end
    total /= count
end

function monte_carlo2(f::Function, D::Int; batch = 50, Z = 1.0)

    means = Real[]

    while true
        points = Real[]
        for i in 1:batch
            x = ntuple(i -> rand(), D)
            push!(points, f(x))
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

function mcerr2(f::Function, g::Function, D::Int; batch = 50, Z = 1.0)
    return sqrt(monte_carlo2(x->(f(x)-g(x))^2, D; batch = batch, Z = Z))
end