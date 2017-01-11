using EllipticalScenGen: MonotonicRiskRegion, SurvivorApproximator
using Distributions
using DataStructures
using StatsFuns
using Base: Test

N = 10000
α = 0.99
β = 0.95
dist = MvNormal(eye(2))
sample = rand(dist, N)
surv = SurvivorApproximator(sample, α)
Ω = MonotonicRiskRegion(surv, 0.95)

count = 0
for i in 1:1000
    x = rand(dist)
    if x ∈ Ω count+=1 end
end
approx_prob = count/1000
