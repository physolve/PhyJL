using Plots, Polynomials
xs = range(0, 10, length = 10) 
ys = @.exp(-xs) 
f = fit(xs, ys) # degree = length(xs) - 1 
f2 = fit(xs, ys, 2) # degree = 2 
print(f2)
scatter(xs, ys, markerstrokewidth = 0, label = "Data") 
plot!(f, extrema(xs)..., label = "Fit") 
plot!(f2, extrema(xs)..., label = "Quadratic Fit")