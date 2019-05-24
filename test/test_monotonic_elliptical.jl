using TailRiskScenGen
using Cones
using LinearAlgebra: I
using Distributions

d = 30
β = 0.95
num_scen = 5000
A = rand(d,d)/sqrt(d)
dist = MvNormal(0.1*rand(d)-0.05, A'A)

K = PolyhedralCone(Array{Float64}(I, d, d))
Ω_gen = EllipticalRiskRegion(dist, K, β)
Ω_mono = MonotonicEllipticalRiskRegion(dist, K, β)
Ω_mono.max_frontier_points = d*5

rights = 0
wrongs = 0
for i in 1:200
    x = rand(dist)
    mon = x ∈ Ω_mono
    gen = x ∈ Ω_gen
    if mon == gen
        rights+=1
    else
        wrongs+=1
    end
    
    if i % 10 == 1
        println("i = ", i, ", num risk frontier points: ", length(Ω_mono.risk_frontier),
              ", num nonrisk frontier points: ", length(Ω_mono.nonrisk_frontier))
    end
end

println("Correct decisions: ", rights)
println("Wrong decisions: ", wrongs)


x = rand(dist)

x ∈ Ω_gen
x ∈ Ω_mono
y = rand(dist)
y ∈ Ω_gen
y ∈ Ω_mono

Ω_mono.risk_frontier
Ω_mono.nonrisk_frontier
Ω=Ω_mono

print("EllipticalRiskRegion...")
tic()
TailRiskScenGen._aggregation_sampling(dist, Ω_gen, num_scen)
toc()

print("MonotonicEllipticalRiskRegion...")
tic()
TailRiskScenGen._aggregation_sampling(dist, Ω_mono, num_scen)
toc()

println("Num non-risk frontier points: ", length(Ω_mono.risk_frontier))
println("Num risk frontier points: ", length(Ω_mono.nonrisk_frontier))
