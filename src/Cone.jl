import Base.length

abstract Cone

@doc """
# Description
Type representing a finitely generated cone. This is the positive
hull of a finite collection of vectors a₁, … , aₙ:

{ ∑ λᵢ xᵢ : λ ≥ 0 }

# Arguments
`A::Matrix{Float64}`: Matrix of cone generators, where each column corresponds to a generator
""" ->
type FiniteCone <: Cone
    A::Matrix{Float64}    # Matrix where each column is a cone generator
    AtA::Matrix{Float64}  # Cross-product of cone generator matrix
    num_gen::Int64
    FiniteCone{T<:Real}(A::Matrix{T}) = new(float(A), float(A'A), size(A,2))
end

function project(cone::FiniteCone, p::Vector{Float64})
    w,z = lcp_solve(cone.AtA, -cone.A'p)
    cone.A * z
end

# Dimension of ambient space
length(cone::FiniteCone) = size(cone.A, 1)

@doc """
# Description
Type representing a Polyhedra cone. This is the
intersection of a finite collection of half-spaces

{x : aᵢ x ≥ 0 ∀ i} = { x: Ax ≥ 0 }

# Arguments
`A::Matrix{Float64}`: Constraint matrix of polyhedral cone
""" ->
type PolyhedralCone{T<:Real} <: Cone
    A::Matrix{T}
    n::Int   # Dimension
    m::Int   # Number of constraints
    PolyhedralCone(A::Matrix{T}) = new(A, size(A,2), size(A,1))
end

PolyhedralCone{T<:Real}(A::Matrix{T}) = PolyhedralCone{T}(A)
function PolyhedralCone{T<:Integer}(A::Matrix{T}, check=false)
    if check
        A = remove_redundant_constraints(A)
    end
    PolyhedralCone{T}(A)
end

function project(cone::PolyhedralCone, p::Vector{Float64})
    res = quadprog(-2.0*p, 2.0*eye(cone.n), cone.A,  '>', zeros(cone.m), fill(-Inf,cone.n), fill(Inf, cone.n),
                   GurobiSolver(OutputFlag=0))
    return res.sol
end

# Dimension of ambient space
length(cone::PolyhedralCone) = size(cone.A, 2)
