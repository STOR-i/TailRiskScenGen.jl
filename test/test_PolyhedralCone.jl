# Test whether PolyhedralCone projects points onto same locations as equivalent FiniteCones

using Distributions

srand(1)

print("\tTest 1 - dim == 2...")
Z = Normal()
μ = rand(Z, 2)
A = rand(Z, (2, 2))
Σ = A'A
dist = MvNormal(μ, Σ)

K_f = RiskRegion(dist, FiniteCone(eye(2)), 0.95).K
K_p = RiskRegion(dist, PolyhedralCone(eye(2)), 0.95).K

for i in 1:1000
    x = rand(dist)
    @test_approx_eq_eps project(K_f, x) project(K_p, x) 1e-5
end
