module EllipticalScenGen

VERSION < v"0.4-" && using Docile

@document

using Distributions, SkewDist
using MathProgBase
using Gurobi
using Clustering

import PDMats: dim

export aggregate_scenarios, Cone, FiniteCone, PolyhedralCone, chernikova, chernikova_general, lcp_solve, project, RiskRegion, transformed_size, in_RiskRegion, aggregation_sampling, nonrisk_clustering, dim, cone_from_constraints, quota_cone

include("LCP_julia.jl")
include("chernikova.jl")
include("chernikova_general.jl")
include("Cone.jl")
include("utilities.jl")
include("RiskRegion.jl")
include("nonrisk_clustering.jl")
include("BruteForceRiskRegion.jl")
include("ScenarioReduction.jl")

end # module
