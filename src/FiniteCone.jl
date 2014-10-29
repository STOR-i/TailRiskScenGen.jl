type FiniteCone
    A::Matrix{Float64}    # Matrix where each column is a cone generator
    AtA::Matrix{Float64}  # Cross-product of cone generator matrix
    dim::Int64
    num_gen::Int64
    FiniteCone(A::Matrix{Float64}) = new(A, A'A, size(A,1), size(A,2))
end

function project(cone::FiniteCone, p::Vector{Float64})
    w,z = lcp_solve(cone.AtA, -cone.A'p)
    cone.A * z
end
