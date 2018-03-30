module TailRiskScenGen

using Distributions
# using SkewDist
using Cones
using Gurobi
using Clustering
using DataStructures
using StatsFuns

import Base: ∈, length
import PDMats: dim

export aggregate_scenarios, Cone, FiniteCone, PolyhedralCone, project, AbstractRiskRegion, EllipticalRiskRegion, MonotonicEllipticalRiskRegion, transformed_size, in_RiskRegion, aggregation_sampling, nonrisk_clustering, dim, cone_from_constraints, quota_cone, MonotonicRiskRegion, SurvivorApproximator

const VecF64 = Vector{Float64}
const MatF64 = Matrix{Float64}

include("utilities.jl")
include("RiskRegion.jl")
include("monotonic_utils.jl")
include("monotonic_risk_regions.jl")
include("nonrisk_clustering.jl")
include("BruteForceRiskRegion.jl")
include("ScenarioReduction.jl")

end # module
