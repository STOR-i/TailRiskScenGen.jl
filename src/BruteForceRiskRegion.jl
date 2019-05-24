mutable struct BruteForceRiskRegion
    μ::Vector{Float64}
    Σ::Matrix{Float64}
    lattice::Array{Vector{Float64}}
    α::Float64
    function BruteForceRiskRegion(μ::Vector{Float64}, Σ::Matrix{Float64}, K::FiniteCone, α::Float64, lattice_width::Int64 = 100)
        checkRiskRegionArgs(μ, Σ, K, α)
        cone_coords = simplex(lattice_width, K.num_gen)
        num_points = length(cone_coords)
        lattice = Array{Vector{Float64}}(undef, num_points)
        for (i,p) in enumerate(cone_coords)
            lattice[i] = -K.A*p
        end
        return new(μ, Σ, lattice, α)
    end
end

function BruteForceRiskRegion(dist::MvNormal, K::FiniteCone, β::Float64, lattice_width::Int64 = 100)
    BruteForceRiskRegion(dist.μ, dist.Σ.mat, K, quantile(Normal(), β), lattice_width)
end

function BruteForceRiskRegion(dist::MvTDist, K::FiniteCone, β::Float64, lattice_width::Int64 = 100)
    BruteForceRiskRegion(dist.μ, dist.Σ.mat, K, quantile(TDist(dist.df), β), lattice_width)
end


function ∈(y::Vector{Float64}, Ω::BruteForceRiskRegion)
    for x in Ω.lattice
        if dot(x,y-Ω.μ) >= sqrt(dot(x,Ω.Σ*x)) * Ω.α
            return true
        end
    end
    false
end    

##
# Returns a matrix whose columns are points of a simplex
# Arguments:
# N - the sum of the elements of each point in the simplex
# n -- the dimension of the space containing the simplex
##
function simplex(N::Int64, n::Int64)
    simplices = Array{Vector{Int64}}(undef, 0)
    if n > 1
        for i in 0:N
            for p in simplex(N-i, n-1)
                push!(p, i)
                push!(simplices, p)
            end
        end
        return simplices
    end
    if n == 1
        points = Array{Vector{Int64}}(undef, 1)
        points[1] = [N]
        return points
    end
end

dim(Ω::BruteForceRiskRegion) = length(Ω.μ)
