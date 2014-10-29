module EllipticalScenGen

export aggregate_scenarios, FiniteCone, lcp_solve, project

include("LCP_julia.jl")
include("FiniteCone.jl")
include("ScenarioReduction.jl")

end # module
