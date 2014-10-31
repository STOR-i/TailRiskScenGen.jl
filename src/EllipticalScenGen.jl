module EllipticalScenGen

export aggregate_scenarios, FiniteCone, lcp_solve, project, RiskRegion, transformed_size, in_RiskRegion, aggregation_sampling

include("LCP_julia.jl")
include("FiniteCone.jl")
include("RiskRegion.jl")
include("ScenarioReduction.jl")

end # module
