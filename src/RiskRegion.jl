function checkRiskRegionArgs(μ::Vector{Float64}, Σ::Matrix{Float64},
                             K::FiniteCone, α::Float64)
        if length(μ) != size(Σ, 1) ArgumentError("μ and Σ must have consistent dimensions") end
        if !isposdef(Σ) ArgumentError("Σ must be a positive definite matrix") end
        if length(K) != length(μ) ArgumentError("FiniteCone must be in same dimension as mean vector") end
        if α <= 0 ArgumentError("α must be strictly positive") end
end
    
type RiskRegion
    μ::Vector{Float64}
    K::Cone
    inv_P::Matrix{Float64}
    α::Float64
    function RiskRegion(μ::Vector{Float64}, Σ::Matrix{Float64},
                        K::FiniteCone, α::Float64)
        checkRiskRegionArgs(μ, Σ, K, α)
        P = chol(Σ)
        inv_P = inv(P)
        new_K = FiniteCone(P*K.A)
        new(μ, new_K, inv_P, α)
    end

    function RiskRegion(μ::Vector{Float64}, Σ::Matrix{Float64},
                        K::PolyhedralCone, α::Float64)
        checkRiskRegionArgs(μ, Σ, K, α)
        P = chol(Σ)
        inv_P = inv(P)
        new_K = PolyhedralCone(K.A * inv_P)
        new(μ, new_K, inv_P, α)
    end
end

RiskRegion(dist::MvNormal, K::FiniteCone, β::Float64) = RiskRegion(dist.μ, dist.Σ.mat, K, quantile(Normal(), β))
RiskRegion(dist::MvTDist, K::FiniteCone, β::Float64) = RiskRegion(dist.μ, dist.Σ.mat, K, quantile(TDist(dist.df), β))

function transformed_size(Ω::RiskRegion, x::Vector{Float64})
    x_trans = x - Ω.μ
    y = Ω.inv_P'x_trans
    norm(project(Ω.K, y))
end

function ∈(x::Vector{Float64}, Ω::RiskRegion)
    transformed_size(Ω::RiskRegion, x::Vector{Float64}) >= Ω.α
end

dim(Ω::RiskRegion) = length(Ω.μ)
