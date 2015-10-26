module EllipticalScenGen

VERSION < v"0.4-" && using Docile

using Distributions, SkewDist
using cones
using Gurobi
using Clustering

import PDMats: dim

export aggregate_scenarios, Cone, FiniteCone, PolyhedralCone, project, RiskRegion, transformed_size, in_RiskRegion, aggregation_sampling, nonrisk_clustering, dim, cone_from_constraints, quota_cone

include("utilities.jl")
include("RiskRegion.jl")
include("nonrisk_clustering.jl")
include("BruteForceRiskRegion.jl")
include("ScenarioReduction.jl")

end # module
