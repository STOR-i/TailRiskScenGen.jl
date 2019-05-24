abstract type AbstractRiskRegion end

@deprecate RiskRegion EllipticalRiskRegion

mutable struct EllipticalRiskRegion <: AbstractRiskRegion
    μ::Vector{Float64}
    K::Cone
    inv_P::Matrix{Float64}
    α::Float64
    function EllipticalRiskRegion(μ::Vector{Float64}, Σ::Matrix{Float64},
                        K::FiniteCone, α::Float64)
        checkRiskRegionArgs(μ, Σ, K, α)
        P = LinearAlgebra.cholesky(Σ).U
        inv_P = inv(P)
        new_K = FiniteCone(P*(-K.A))
        new(μ, new_K, inv_P, α)
    end

    function EllipticalRiskRegion(μ::Vector{Float64}, Σ::Matrix{Float64},
                        K::PolyhedralCone, α::Float64)
        checkRiskRegionArgs(μ, Σ, K, α)
        P = LinearAlgebra.cholesky(Σ).U
        inv_P = inv(P)
        new_K = PolyhedralCone{Float64}((-K.A) * inv_P) # quickfix
        new(μ, new_K, inv_P, α)
    end
end

function EllipticalRiskRegion(μ::Distributions.ZeroVector, Σ::Matrix{Float64}, K::Cone, α::Float64)
    EllipticalRiskRegion(zeros(μ.len), Σ, K, α)
end
EllipticalRiskRegion(dist::AbstractMvNormal, K::Cone, β::Float64) = EllipticalRiskRegion(dist.μ, dist.Σ.mat, K, quantile(Normal(), β))
EllipticalRiskRegion(dist::Distributions.AbstractMvTDist, K::Cone, β::Float64) = EllipticalRiskRegion(dist.μ, dist.Σ.mat, K, quantile(TDist(dist.df), β))
# EllipticalRiskRegion(dist::MvSkewTDist, K::Cone, β::Float64) = EllipticalRiskRegion(dist.ξ, dist.Ω.mat, K, quantile(TDist(dist.df), β))

function transformed_size(Ω::EllipticalRiskRegion, x::Vector{Float64})
    x_trans = x - Ω.μ
    y = Ω.inv_P'x_trans
    norm(project(Ω.K, y))
end

function ∈(x::Vector{Float64}, Ω::EllipticalRiskRegion)
    transformed_size(Ω::EllipticalRiskRegion, x::Vector{Float64}) >= Ω.α
end

length(Ω::EllipticalRiskRegion) = length(Ω.μ)

mutable struct MonotonicEllipticalRiskRegion <: AbstractRiskRegion
    # Note that cost function for portfolio selection is negative
    # monotonic, so roles of risk and non-risk frontiers are reversed
    Ω::EllipticalRiskRegion
    nonrisk_frontier::LinkedList{VecF64}
    risk_frontier::LinkedList{VecF64}
    max_frontier_points::Int
    function MonotonicEllipticalRiskRegion(args...)
        Ω=EllipticalRiskRegion(args...)
        new(Ω, nil(VecF64), nil(VecF64), length(Ω)*10)
    end
end

function ∈(x::VecF64, Ω::MonotonicEllipticalRiskRegion)
    # Note that cost function for portfolio selection is negative
    # monotonic, so roles of risk and non-risk frontiers are reversed
    if above_risk_frontier(x, Ω.risk_frontier) return false
    elseif below_nonrisk_frontier(x, Ω.nonrisk_frontier) return true
    else
        if x ∈ Ω.Ω
            if length(Ω.nonrisk_frontier) < Ω.max_frontier_points
                add_to_nonrisk_frontier!(Ω, x)
            end
            return true
        else
            if length(Ω.risk_frontier) < Ω.max_frontier_points
                add_to_risk_frontier!(Ω, x)
            end
            return false
        end
    end
end
