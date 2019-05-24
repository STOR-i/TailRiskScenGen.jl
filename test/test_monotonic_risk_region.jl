using TailRiskScenGen: MonotonicRiskRegion, SurvivorApproximator, prob_nonrisk
using Distributions
using DataStructures
using LinearAlgebra: I
using StatsFuns

d = 2
N = 10000
α = 0.99
β = 0.95
dist = MvNormal(Array{Float64}(I, d, d))
sample = rand(dist, N)
surv = SurvivorApproximator(sample, α)
Ω = MonotonicRiskRegion(surv, 0.95)

x = rand(dist)
x ∈ Ω

tic()
count = 0
for i in 1:1000
    x = rand(dist)
    if !(x ∈ Ω) count+=1 end
end
approx_prob = count/1000
toc()

α = 0.99
β_list = [0.95, 0.99]
num_points_nonrisk = 1000
num_points_surv = 10000
dist = MvNormal(eye(d))
prob_nonrisk(dist, β_list, num_points_surv, num_points_nonrisk, α)
