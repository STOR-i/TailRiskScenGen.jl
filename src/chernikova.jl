# Tests whether two vectors are collinear within a certain tolerance
function collinear(a::Array{Float64}, b::Array{Float64}, tol::Float64 = 1e-6)
    length(a) == length(b) || error("Input arrays must have the same dimenisions")
    size_a, size_b = norm(a), norm(b)
    if size_a == 0.0 || size_b == 0.0; error("Input vectors must be non-zero"); end
    a_n, b_n = a/norm(a), b/norm(b)
    if norm(a_n - b_n) < tol || norm(a_n + b_n) < tol; return true end;
    return false
end

type ChernMat{T<:Real}
    B::Matrix{T}
    m::Int64     # Number of constraints
    n::Int64     # Dimension
    rows::Int    # Current number of rows
end

function ChernMat{T<:Real}(A::Matrix{T})
    m, n = size(A)
    ChernMat([eye(T,n) A'], m, n, n)
end

function pprint(mat::ChernMat{Float64})
    @printf("Rows: %7d\n", mat.rows)
    for i in 1:mat.rows
        for j in 1:mat.n
            @printf("%7.3f", mat.B[i,j])
        end
        print(" |")
        for j in mat.n+1:mat.n+mat.m
            @printf("%7.3f", mat.B[i,j])
        end
        print("\n")
    end
end

function pprint(mat::ChernMat{Int})
    @printf("Rows: %7d\n", mat.rows)
    for i in 1:mat.rows
        for j in 1:mat.n
            @printf("%6d", mat.B[i,j])
        end
        print(" |")
        for j in mat.n+1:mat.n+mat.m
            @printf("%6d", mat.B[i,j])
        end
        print("\n")
    end
end


# Main algorithm
# Calculates the finite generators of the cone formed
# by the intersection of a polyhedral cone and the positive
# quadrant; that is the cone formed by the positive solutions
# of the following problem:
#    Ax >= 0
#
#
# This is solved by Chernikova's algorithm:
# "Algorithm for finding a general formula for the non-negative
# solutions of linear inequalities" by N.V. Chernikova (1964)
#
# Returns a matrix whose columns are the required cone generators
function chernikova{T<:Real}(A::Matrix{T}, verbosity::Int = 0)
    mat = ChernMat(A)
    counter = 1
    while true
        if verbosity > 0
            println("Iteration ", counter)
            println("---------------")
            pprint(mat)
        end
        (col, p) = leading_col(mat)
        if col == 0
            if verbosity > 0; println("Generators have been found!"); end;
            return LHS(mat)'
        elseif p == mat.rows
            if verbosity > 0; println("Zero is the unique solution"); end;
            return zeros(n)
        end
        
        # Inspect matrix to find the maximum size of the next
        if verbosity > 0; println("Using column $(col) which has $(p) negatives"); end;
        pos, nil, neg = get_pos_neg_indices(mat, col)
        max_rows = length(pos) + length(nil) + length(pos) * length(neg)
        new_B = Array(T, max_rows, mat.m + mat.n)


        # Handle case of two rows separately
        if mat.rows == 2
            if verbosity > 0; println("Handling two row case..."); end;
            new_B[1,:] = mat.B[pos[1],:]
            new_row = combine(mat, col, pos[1], neg[1])
            if verbosity > 1; println("Candidate row: ", new_row); end;
            if !collinear(mat.B[pos[1],:], new_row)
                new_B[2,:] = new_row
                mat.rows = 2
                mat.B = new_B
            else
                mat.rows = 1
                mat.B = new_B[1,:]
            end
            continue
        end
                
        # Add rows with non-negative intersections
        # with leading column to next matrix
        for (k,i) in enumerate([pos,nil]); new_B[k,:] = mat.B[i,:]; end;
        new_rows = length(pos) + length(nil)

        # Go through pairs of row with positive and negative
        # coefficients in leading column and combine if necessary
        for i1 in pos, i2 in neg
            if combinable(mat, i1, i2, verbosity - 1)
                new_rows += 1
                new_B[new_rows, :] = combine(mat, col, i1, i2)
            end
        end
        mat.B = new_B[1:new_rows,:]
        mat.rows = new_rows
        counter += 1

        if verbosity > 0
            @printf("Used %d rows with non-negative intersections with leading column\n", length(pos) + length(nil))
            @printf("Added %d positive combinations of rows out of a possible %d\n",
                    new_rows - length(pos) - length(nil), length(pos) * length(neg))
            println("------------------------------------------------")
        end
    end
end

function leading_col{T<:Real}(mat::ChernMat{T})
    p = mat.rows
    I = 0
    for i in mat.n+1:mat.n+mat.m
        r = sum(mat.B[:,i] .< 0.0)
        if r == mat.rows
            return (i, r)
        elseif (r > 0) && (r < p)
           I, p = i, r
        end
    end
    return (I, p)
end

LHS(mat::ChernMat) = mat.B[:, 1:mat.n]

function get_pos_neg_indices{T<:Real}(mat::ChernMat{T}, col::Int64)
    pos, nil, neg = Int[],Int[],Int[]
    for k in 1:mat.rows
        if mat.B[k,col] > 0.0; push!(pos, k)
        elseif mat.B[k,col] < 0.0; push!(neg, k)
        else push!(nil, k) end
    end
    pos, nil, neg
end
    
function combinable{T<:Real}(mat::ChernMat{T}, i1::Int64, i2::Int64, verbosity::Int = 0)
    # Find all columns for which the rows have a common zero
    # LHS
    zero_cols = Int[]
    for j in 1:mat.n
        if (mat.B[i1,j] == 0.0) && (mat.B[i2,j] == 0.0)
            push!(zero_cols, j)
        end
    end

    # RHS
    for j in mat.n+1:mat.n + mat.m
        if all(mat.B[:,j] .>= 0.0)
            if (mat.B[i1,j] == 0.0) && (mat.B[i2,j] == 0)
                push!(zero_cols, j)
            end
        end
    end
    if verbosity > 0;
        println("Rows $(i1) and $(i2) have ", length(zero_cols), " permissible common zeros at columns ", zero_cols);
    end;
    if length(zero_cols) == 0 return false end
    for i in 1:mat.rows
        if i == i1 || i == i2; continue; end
        
        if verbosity > 1; println("\tRow $i - Values of zero cols: ",  mat.B[i,zero_cols]); end;
        if vec(mat.B[i,zero_cols]) == zeros(length(zero_cols))
            return false
        end
    end
    return true
end


function combine{T<:Real}(mat::ChernMat{T}, col::Int, i1::Int, i2::Int)
    a, b = 1.0, -mat.B[i1,col]/mat.B[i2,col]
    return a*mat.B[i1,:] + b*mat.B[i2,:]
end

function combine(mat::ChernMat{Int}, col::Int, i1::Int, i2::Int)
    t = lcm(mat.B[i1, col], mat.B[i2, col])
    return (t/mat.B[i1,col]) * mat.B[i1, :] - (t/mat.B[i2,col]) * mat.B[i2, :]
end

# Find the conical hull defined by the following constraints:
#   1ᵀx = c
#   Ax ≤ b
#   x ≥ 0
#  Returns FiniteCone object
function cone_from_constraints(A::Matrix{Float64}, b::Vector{Float64}, c::Float64)
    m, n = size(A)
    m == length(b) || error("A and b must have consistent dimensions")
    c > 0 || error("c must be strictly positive")
    A_poly = Array(Float64, m+1, n)
    A_poly[1,:] = ones(n)
    for i in 1:m
        A_poly[i+1,:] = (b[i]/c) * ones(n)' - A[i,:]
    end
    
    return FiniteCone(chernikova(A_poly))
end

# Finds the conical hull of of the region define
# by the following constraints:
#
# 1ᵀ x = 1
# x ≤ u
# x ≥ 0
# Returns FiniteCone object
function quota_cone(u::Vector{Float64})
    n = length(u)
    return cone_from_constraints(eye(n), u, 1.0)
end

