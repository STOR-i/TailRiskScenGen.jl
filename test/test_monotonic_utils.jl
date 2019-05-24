using TailRiskScenGen: vec_isless, SurvivorApproximator, below_nonrisk_frontier, above_risk_frontier, decapitate!
using Distributions
using LinearAlgebra: I
using DataStructures

A = [[1.0, 2.0, 3.0] [2.0, 0.0, 6.0] [3.0, 8.0, 9.0]]
@test vec_isless([1.0, 2.0, 3.0], [2.0, 3.0, 4.0])
@test !vec_isless([1.0, 1.0], [1.0, 2.0])
@test !vec_isless(view(A,:,1), view(A,:,2))
@test vec_isless(view(A,:,1), view(A,:,3))
@test !vec_isless(view(A,:,3), [10.0, 9.0, 7.0])

d = 3
dist = MvNormal(Array{Float64}(I, d, d))
sample = rand(dist, 10000)
α = 0.99
surv = SurvivorApproximator(sample, α)
approx_prob, err = surv(zeros(d))
@test 0.5^d ≈ approx_prob atol=err

frontier=list([2.0, 0.0, 1.0], [1.0, 1.0, 1.0], [1.0, 0.0, 2.0])
@test above_risk_frontier([3.0, 3.0, 3.0], frontier)
@test above_risk_frontier([1.5,1.5,1.5], frontier)
@test !above_risk_frontier([3.0, 0.0, 3.0], frontier)
@test !above_risk_frontier([-2.0, -1.0, -3.0], frontier)
@test below_nonrisk_frontier([0.0, 0.0, 0.0], frontier)
@test below_nonrisk_frontier([0.5, -0.1, 1.5], frontier)
@test !below_nonrisk_frontier([1.0, 1.0, 1.0], frontier)
@test !below_nonrisk_frontier([5.0, 5.0, 1.0], frontier)

l = list(1,2,3,4)
l2 = tail(l)
decapitate!(l2)
@test l == list(1,3,4)
@test l2 == list(3,4)

