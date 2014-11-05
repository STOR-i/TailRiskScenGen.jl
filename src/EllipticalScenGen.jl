module EllipticalScenGen

using Distributions

export aggregate_scenarios, FiniteCone, lcp_solve, project, RiskRegion, transformed_size, in_RiskRegion, aggregation_sampling, nonrisk_clustering

include("LCP_julia.jl")
include("FiniteCone.jl")
include("RiskRegion.jl")
include("Clustering.jl")
include("BruteForceRiskRegion.jl")
include("ScenarioReduction.jl")

end # module
