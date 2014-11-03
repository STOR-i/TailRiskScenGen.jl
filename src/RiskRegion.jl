type RiskRegion
    μ::Vector{Float64}
    Σ::Matrix{Float64}
    K::FiniteCone
    inv_P::Matrix{Float64}
    α::Float64
    function RiskRegion(μ::Vector{Float64}, Σ::Matrix{Float64},
                        K::FiniteCone, α::Float64)
        if length(μ) != size(Σ, 1) ArgumentError("μ and Σ must have consistent dimensions") end
        if !isposdef(Σ) ArgumentError("Σ must be a positive definite matrix") end
        if K.dim != length(μ) ArgumentError("FiniteCone must have same dimension as mean vector") end
        if α <= 0 ArgumentError("α must be strictly positive") end
        P = chol(Σ)
        inv_P = inv(P)
        new_K = FiniteCone(P*K.A)
        new(μ, Σ, new_K, inv_P, α)
    end
end

function RiskRegion(dist::MvNormal, K::FiniteCone, β::Float64)
    RiskRegion(dist.μ, dist.Σ.mat, K, quantile(Normal(), β))
end

function transformed_size(Ω::RiskRegion, x::Vector{Float64})
    x_trans = x - Ω.μ
    y = Ω.inv_P'x_trans
    norm(project(Ω.K, y))
end

function in_RiskRegion(Ω::RiskRegion, x::Vector{Float64})
    transformed_size(Ω::RiskRegion, x::Vector{Float64}) >= Ω.α
end

