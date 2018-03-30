using TailRiskScenGen
using Base: Test
using Distributions

function valid_scenario_set(scenarios::Matrix{Float64}, probs::Vector{Float64})
    @test size(scenarios, 2) == length(probs)
    @test sum(probs) ≈ 1.0
    for p in probs
        @test p >= 0
    end
end

d = 10
num_scen = 1000
μ = rand(Normal(), d)
A = rand(d,d)./d
Σ = A'A
K = FiniteCone(-eye(d))
dist = MvNormal(μ, Σ)
Ω = EllipticalRiskRegion(dist, K, 0.95)
scenarios = rand(dist, num_scen)

new_scen, new_prob = aggregate_scenarios(scenarios, Ω)
valid_scenario_set(new_scen, new_prob)

agg_sample_scen, agg_sample_prob = aggregation_sampling(dist, Ω, num_scen)
valid_scenario_set(agg_sample_scen, agg_sample_prob)

cluster_sample_scen, cluster_sample_prob = nonrisk_clustering(dist, Ω, 100,10)
valid_scenario_set(cluster_sample_scen, cluster_sample_prob)

cluster_scen, cluster_probs = nonrisk_clustering(scenarios, Ω, 10)
valid_scenario_set(cluster_scen, cluster_probs)

# Monotonic Risk Region
sample = rand(dist, 10000)
α = 0.99
surv = SurvivorApproximator(sample, α)
Ω = MonotonicRiskRegion(surv, 0.95)
new_scen, new_prob = aggregate_scenarios(scenarios, Ω)
valid_scenario_set(new_scen, new_prob)

agg_sample_scen, agg_sample_prob = aggregation_sampling(dist, Ω, num_scen)
valid_scenario_set(agg_sample_scen, agg_sample_prob)

cluster_sample_scen, cluster_sample_prob = nonrisk_clustering(dist, Ω, 100, 10)
valid_scenario_set(cluster_sample_scen, cluster_sample_prob)

cluster_scen, cluster_probs = nonrisk_clustering(scenarios, Ω, 10)
valid_scenario_set(cluster_scen, cluster_probs)
