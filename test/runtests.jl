using TailRiskScenGen
using Test

println("Running test_monotonic_utils.jl...")
include("test_monotonic_utils.jl")

println("Running test_RiskRegion.jl...")
include("test_RiskRegion.jl")

println("Running test_BruteForceRiskRegion.jl...")
include("test_BruteForceRiskRegion.jl")

println("Running test_scen_gen.jl...")
include("test_scen_gen.jl")

println("Done tests")

## # Parameters
## mean = [0.1, 0.05]
## cov = [[1.0, -0.5] [-0.5, 1.5]]
## num_scen = 100
## beta = 0.95

## dim = length(mean)
## mv_norm = MvNormal(mean, cov)
## scenarios = rand(mv_norm, num_scen)
## probs = fill(1.0/num_scen, num_scen)
## alpha = quantile(Normal(), 0.95)
## cone_A = eye(dim)

## print("Running aggregate_scenarios with test arguments...")
## agg_scen, agg_prob = aggregate_scenarios(scenarios, probs, mean, cov, cone_A, alpha)
## print(agg_scen)
## print(agg_prob)
## print("\n")
