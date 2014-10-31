using Distributions

function valid_scenario_set(scenarios::Matrix{Float64}, probs::Vector{Float64})
    @test size(scenarios, 2) == length(probs)
    @test_approx_eq sum(probs) 1.0
    for p in probs
        @test (p >= 0)
    end
end

dim = 10
num_scen = 1000
μ = rand(Normal(), dim)
A = rand(Normal(), (dim,dim))
Σ = A'A
K = FiniteCone(eye(dim))
Ω = RiskRegion(μ, Σ, K, 1.96)
scenarios = rand(MvNormal(μ, Σ), num_scen)

new_scen, new_prob = aggregate_scenarios(scenarios, Ω)
valid_scenario_set(new_scen, new_prob)
