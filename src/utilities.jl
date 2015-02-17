# Find the conical hull defined by the following constraints:
#   1ᵀx = c
#   Ax ≤ b
#   x ≥ 0
#  Returns PolyhedralCone object
function cone_from_constraints(A::Matrix{Float64}, b::Vector{Float64}, c::Float64)
    m, n = size(A)
    m == length(b) || error("A and b must have consistent dimensions")
    c > 0 || error("c must be strictly positive")
    A_poly = Array(Float64, m+1, n)
    A_poly[1,:] = ones(n)
    for i in 1:m
        A_poly[i+1,:] = (b[i]/c) * ones(n)' - A[i,:]
    end
    
    return PolyhedralCone(A_poly)
end

# Finds the conical hull of of the region define
# by the following constraints:
#
# 1ᵀ x = 1
# x ≤ u
# x ≥ 0
# Returns PolyhedralCone object
function quota_cone(u::Vector{Float64})
    n = length(u)
    return cone_from_constraints(eye(n), u, 1.0)
end

function checkRiskRegionArgs(μ::Vector{Float64}, Σ::Matrix{Float64},
                             K::Cone, α::Float64)
        if length(μ) != size(Σ, 1) ArgumentError("μ and Σ must have consistent dimensions") end
        if !isposdef(Σ) ArgumentError("Σ must be a positive definite matrix") end
        if length(K) != length(μ) ArgumentError("Cone must be in same dimension as mean vector") end
        if α <= 0 ArgumentError("α must be strictly positive") end
end
    
