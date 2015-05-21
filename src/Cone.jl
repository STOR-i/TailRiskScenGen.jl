import Base.length

abstract Cone

type FiniteCone <: Cone
    A::Matrix{Float64}    # Matrix where each column is a cone generator
    AtA::Matrix{Float64}  # Cross-product of cone generator matrix
    num_gen::Int64
    FiniteCone(A::Matrix{Float64}) = new(A, A'A, size(A,2))
end

function project(cone::FiniteCone, p::Vector{Float64})
    w,z = lcp_solve(cone.AtA, -cone.A'p)
    cone.A * z
end

# Dimension of ambient space
length(cone::FiniteCone) = size(cone.A, 1)

type PolyhedralCone <: Cone
    A::Matrix{Float64}
    n::Int   # Dimension
    m::Int   # Number of constraints
    PolyhedralCone(A::Matrix{Float64}) = new(A, size(A,2), size(A,1))
end

function project(cone::PolyhedralCone, p::Vector{Float64})
    res = quadprog(-2.0*p, 2.0*eye(cone.n), cone.A,  '>', zeros(cone.m), fill(-Inf,cone.n), fill(Inf, cone.n),
                   GurobiSolver(OutputFlag=0))
    return res.sol
end

# Dimension of ambient space
length(cone::PolyhedralCone) = size(cone.A, 2)
