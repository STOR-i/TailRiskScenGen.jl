function collinear(a::Array{Int}, b::Array{Int})
    length(a) == length(b) || error("Input arrays must have the same dimenisions")
    a_n, b_n = a/gcd(a), b/gcd(b)
    if a_n==b_n || a_n==-b_n; return true; end
    return false
end

# Normalises the columns of a matrix
function norm_cols(A::Matrix{Int})
    norm_A = copy(A)
    for i in 1:size(A,2)
        if (g = gcd(A[:,i])) > 1; norm_A[:, i] = broadcast(div, norm_A[:, i], g); end
    end
    return norm_A
end

type ChernMat{T<:Integer}
    B::Matrix{T}
    m::Int64     # Number of constraints
    n::Int64     # Dimension
    rows::Int    # Current number of rows
end

function ChernMat{T<:Integer}(A::Matrix{T})
    m, n = size(A)
    ChernMat([eye(T,n) A'], m, n, n)
end

function remove_row{T<:Integer}(A::Matrix{T}, k::Int)
    m, n = size(A)
    B = Array(T, m - 1, n)
    i = 1
    for j in 1:m
        if j != k
            B[i,:] = A[j,:]
            i+=1
        end
    end
    return B
end

function add_row{T<:Integer}(A::Matrix{T}, a::Matrix{T})
    # Here a is a row vector
    m, n = size(A)
    B = Array(T, m+1, n)
    B[1:m, 1:n] = A
    B[m+1, :] = a
    return B
end

function pprint(mat::ChernMat)
    @printf("# Rays: %7d\n", mat.rows)
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

function chernikova{T<:Integer}(mat::ChernMat{T}, verbosity::Int = 0)
    m, n = mat.m, mat.n
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
            return norm_cols(LHS(mat)')
        elseif p == mat.rows
            if verbosity > 0; println("Zero is the unique solution"); end;
            return zeros(n)
        end
        
        # Inspect matrix to find the maximum size of the next
        #if verbosity > 0; println("Using column $(col) which has $(p) negatives"); end;
        pos, nil, neg = get_pos_neg_indices(mat, col)
        
        if verbosity > 0
            println("Adding constraint $(col) which is currently violated by $(length(neg)) of $(mat.rows) extremal rays")
            println("New basis will consist of a maximum of $(mat.rows) - $(length(neg)) + $(length(pos)) * $(length(neg)) = $(mat.rows - length(neg) + length(pos)*length(neg)) extremal rays")
        end;
        
        max_rows = length(pos) + length(nil) + length(pos) * length(neg)
        new_B = Array(T, max_rows, mat.m + mat.n)


        # Handle case of two rows separately
        if mat.rows == 2
            if verbosity > 0; println("Handling two row case..."); end;
            new_B[1,:] = mat.B[pos[1],:]
            new_row = combine(mat.B, col, pos[1], neg[1])
            if verbosity > 1; println("Candidate row: ", new_row); end;
            if !collinear(vec(mat.B[pos[1],:]), vec(new_row))
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
                new_B[new_rows, :] = combine(mat.B, col, i1, i2)
            end
        end
        mat.B = new_B[1:new_rows,:]
        mat.rows = new_rows
        counter += 1

        if verbosity > 0
            @printf("Added %d / %d positive combinations of extremal rays\n",
                    new_rows - length(pos) - length(nil), length(pos) * length(neg))
            println("------------------------------------------------")
        end
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
function chernikova{T<:Integer}(A::Matrix{T}, verbosity::Int = 0)
    mat = ChernMat(A)
    chernikova(mat, verbosity)
end

function leading_col{T<:Integer}(mat::ChernMat{T})
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

function get_pos_neg_indices{T<:Integer}(mat::ChernMat{T}, col::Int64)
    pos, nil, neg = Int[],Int[],Int[]
    for k in 1:mat.rows
        if mat.B[k,col] > 0.0; push!(pos, k)
        elseif mat.B[k,col] < 0.0; push!(neg, k)
        else push!(nil, k) end
    end
    pos, nil, neg
end
    
function combinable{T<:Integer}(mat::ChernMat{T}, i1::Int64, i2::Int64, verbosity::Int = 0)
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
    
    if length(zero_cols) == 0
        if verbosity > 0; println("Rays $(i1) and $(i2) saturate no common constraint - they do not combine to form a new extremal ray") end
        return false
    end
        
    if verbosity > 0;
        println("Rays $(i1) and $(i2) have ", length(zero_cols), " are both saturated by constraints ", zero_cols);
    end;
    
    for i in 1:mat.rows
        if i == i1 || i == i2; continue; end
        
        if verbosity > 1; println("\tRay $i - values on saturated constraints: ",  mat.B[i,zero_cols]); end;
        if vec(mat.B[i,zero_cols]) == zeros(length(zero_cols))
            if verbosity > 0; println("\tRay $i saturates all saturated constraints of these rays so do not combine") end;
            return false
        end
    end
    if verbosity > 0;
        println("\tNo other ray saturates the common saturated constraints of these rays so combine to form a new extreme ray")
    end;
    return true
end

function combine{T<:Integer}(mat::Matrix{T}, col::Int, i1::Int, i2::Int)
    a, b = 1.0, -mat[i1,col]/mat[i2,col]
    return a*mat[i1,:] + b*mat[i2,:]
end

function combine{T<:Integer}(mat::Matrix{T}, col::Int, i1::Int, i2::Int)
    t = gcd(mat[i1, col], mat[i2, col])
    a1, a2 = div(mat[i1, col],t), div(mat[i2, col],t)
    return a1 * mat[i2,:] - a2 * mat[i1,:] 
end

function combine{T<:Integer}(col::Int, vec1::Vector{T}, vec2::Vector{T})
    t = gcd(vec1[col], vec2[col])
    a1, a2 = div(vec1[col],t), div(vec2[col],t)
    return a1 * vec2 - a2 * vec1
end
