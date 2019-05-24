# Tests whether risk regions calculate the correct transformed distance
# Results come from previous Python implementation of Risk regions
using LinearAlgebra: diagm
using DelimitedFiles: readdlm

function test_points(Ω::EllipticalRiskRegion, file_name::String)
    res = readdlm(file_name)
    for i in 1:size(res,1)
        row = res[i,:]
        x, dist = row[1:end-1], row[end]
        @test transformed_size(Ω, x) ≈ dist atol=1e-2
    end
end

print("\tTest 1...")
K₁ = FiniteCone(-Array{Float64}(I, 2, 2))
μ₁ = [0.0, 0.0]
Σ₁ = [[2.0, 0.0] [0.0, 2.0]]
Ω₁ = EllipticalRiskRegion(μ₁, Σ₁, K₁, 1.0)
test_points(Ω₁, "results_1.txt")
print("success!\n")

print("\tTest 2...")
K₂ = FiniteCone(-Array{Float64}(I, 3, 3))
μ₂ = [1.0, -0.5, 3.0]
Σ₂ = diagm(0=>[1.0, 2.0, 4.0])
Ω₂ = EllipticalRiskRegion(μ₂, Σ₂, K₂, 1.0)
test_points(Ω₂, "results_2.txt")
print("success!\n")

print("\tTest 3...")
K₃ = FiniteCone(-Array{Float64}(I, 3, 3))
μ₃ = [1.0, -0.5, 3.0]
Σ₃ = [2.0 0.5 0.0;
      0.5 2.5 0.0;
      0.0 0.0 4.0]
Ω₃ = EllipticalRiskRegion(μ₃, Σ₃, K₃, 1.0)
test_points(Ω₃, "results_3.txt")
print("success!\n")

print("\tTest 5...")
K₅ = FiniteCone(-[[1.0, 1.0, 0.0] [0.0, 1.0, -0.5] [1.0, 0.0, 1.0]])
μ₅ = [1.0, -0.5, 3.0]
Σ₅ = [2.0 0.5 0.0;
      0.5 2.5 0.0;
      0.0 0.0 4.0]
Ω₅ = EllipticalRiskRegion(μ₅, Σ₅, K₅, 1.0)
test_points(Ω₅, "results_5.txt")
print("success!\n")


