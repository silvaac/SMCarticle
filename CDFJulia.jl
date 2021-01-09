using Interpolations
using Plots

function rCDF(x,y,n)
    ru = rand(n)
    a = LinearInterpolation(y, x, extrapolation_bc = Flat())
    xout = a(ru)
    return xout
end

function breakCDF(x,y)
    xp =  x .>= 0
    xn =  x .< 0
    yp =  1 .-y[xp]
    yn =  y[xn]
    return x[xp], yp, x[xn] ,yn
end

function getCDF(x)
    sx = sort(x)
    n  = length(x)
    y  = (1:n)/(n+1)
    return sx,y
end

function mountainPlot(x)
    x,y = getCDF(x)
    xp,yp,xn,yn = breakCDF(x,y)
    plot(xp,yp,color="blue",yaxis=:log);plot!(xn,yn,color="blue",yaxis=:log)
end


function getrnd(r,n)
    x,y = getCDF(r)
    xout = rCDF(x,y,n)
    return xout
end
