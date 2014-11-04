# Tests if RiskRegion and BruteForceRiskRegion define (approximately) the same region by
# testing for membership of both for many randomly generated points.

using Distributions

# Change test inputs as desired
dim = 3
num_points = 1000
lattice_width = 100
β = 0.99

#
srand(1)
Z = Normal()
μ = rand(Z, dim)
A = rand(Z, (dim, dim))
Σ = A'A
## μ = fill(0.0, dim)
## Σ = eye(dim)
dist = MvNormal(μ, Σ)
K = FiniteCone(eye(dim))

Ω_ellipse = RiskRegion(dist, K, β)
Ω_bf = EllipticalScenGen.BruteForceRiskRegion(dist, K, β, lattice_width)

# println("μ = ", μ)
# println("Σ = \n", Σ)

for i in 1:num_points
    y = rand(dist)
    ellipse_res = in_RiskRegion(Ω_ellipse, y)
    bf_res = in_RiskRegion(Ω_bf, y)
    if  ellipse_res != bf_res
        println("\tERROR: failed at point ", i, " = ", y)
        println("\tEllipse results = ", ellipse_res)
        println("\tTransformed size = ", transformed_size(Ω_ellipse, y))
        println("\tΩ.α = ", Ω_ellipse.α)
        println("\tBrute Force results = ", bf_res)
        break
    end
end
