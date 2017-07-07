# Benchmarks for projection time of FiniteCone vs. PolyhedralCone

using BenchmarkLite
using Distributions, TailRiskScenGen

# FiniteCone Benchmark test definition 

type ConeProject <: Proc
    dist::ContinuousMultivariateDistribution
    K::Cone
    function ConeProject(dist::ContinuousMultivariateDistribution, K::Cone)
        length(dist)==length(K) || throw(ArgumentError("Distribution and cone must be defined in same dimension"))
        new(dist, K)
    end
end

function Base.string(proc::ConeProject)
    dim = length(proc.K)
    ty = typeof(proc.K)
    "$(ty): Dim=$(dim)"
end
Base.length(proc,cfg) = cfg
Base.isvalid(proc::ConeProject, cfg) = true
function Base.start(proc::ConeProject, cfg)
    rand(proc.dist, cfg)
end

function Base.run(proc::ConeProject, cfg, s)
    for i in 1:cfg
        project(proc.K, s[:,i])
    end
end
Base.done(proc::ConeProject, cfg, s) = nothing

# Run test
dims = [2, 5, 10, 20]
procs = Proc[ConeProject(MvNormal(eye(40)), FiniteCone(eye(40))), ConeProject(MvNormal(eye(40)), PolyhedralCone(eye(40)))]
cfgs = [100, 1000, 2000]

rtable = run(procs, cfgs)
show(rtable)
