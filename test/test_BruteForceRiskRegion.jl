# Tests if RiskRegion and BruteForceRiskRegion define (approximately) the same region by
# testing for membership of both for many randomly generated points.

using Distributions

function verify_RiskRegion(dist::Sampleable, Ω::RiskRegion, Ω_bf::EllipticalScenGen.BruteForceRiskRegion, n::Int64)
    for i in 1:n
        y = rand(dist)
        ellipse_res = y ∈ Ω
        bf_res = y ∈ Ω_bf
        if  ellipse_res != bf_res
            println("\tERROR: failed at point ", i, " = ", y)
            println("\tEllipse results = ", ellipse_res)
            println("\tTransformed size = ", transformed_size(Ω, y))
            println("\tΩ.α = ", Ω.α)
            println("\tBrute Force results = ", bf_res)
            break
        end
    end
end

srand(1)

# Test RiskRegion for multivariate normal distribution

print("\tTest 1 - MvNormal...")
dim = 3
num_points = 1000
lattice_width = 100
β = 0.99


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

verify_RiskRegion(dist, Ω_ellipse, Ω_bf, num_points)

print("success\n")

# Test RiskRegion for multivariate T distribution

print("\tTest 2 - MvTDist...")
dim = 3
num_points = 1000
lattice_width = 100
β = 0.99

Z = Normal()
μ = rand(Z, dim)
A = rand(Z, (dim, dim))
Σ = A'A
df = 4.0
## μ = fill(0.0, dim)
## Σ = eye(dim)
dist = MvTDist(df, μ, Σ)
K = FiniteCone(eye(dim))
Ω_ellipse = RiskRegion(dist, K, β)
Ω_bf = EllipticalScenGen.BruteForceRiskRegion(dist, K, β, lattice_width)

verify_RiskRegion(dist, Ω_ellipse, Ω_bf, num_points)
print("success\n")
# println("μ = ", μ)
# println("Σ = \n", Σ)

