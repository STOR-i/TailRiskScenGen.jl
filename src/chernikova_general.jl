function redundant_constraint_check{T<:Real}(mat::ChernMat{T}, k::Int, verbosity::Int = 0)
    zeros = find(x->x==0, mat.B[:, k+mat.n])
    red_rows = IntSet()
    for i in 1:k-1
        zs = find(x->x==0, mat.B[:, i+mat.n])
        if issubset(zs, zeros) && length(setdiff(zeros, zs)) > 0
            push!(red_rows, i)
        elseif issubset(zeros, zs) && length(setdiff(zs, zeros))
            push!(red_rows, k)
            break
        end
    end
    return red_rows
end

function chernikova_general{T<:Real}(A::Matrix{T}, verbosity::Int = 0)
    # Initialisation
    m, n = size(A)
    bidray = ChernMat(A)
    uniray = ChernMat(Array(T, 0, m+n), m, n, 0)

    redundant_rows=IntSet()
    
    for k in 1:m
        if verbosity > 0
            println("Iteration ", k)
            println("---------------")
            println("Bidirectional rays")
            pprint(bidray)
            println("Unidirectional rays")
            pprint(uniray)
        end
        # Find if it exists, bidirectional ray which does not saturate k-th constraint
        index_non_zero = -1
        for i in 1:bidray.rows
            if bidray.B[i,k+n] != 0
                index_non_zero = i
                break
            end
        end

        if index_non_zero != -1
            if verbosity > 0
                println("Bidirectional ray $(index_non_zero) is unsaturated by constraint $(k)")
            end
            new_uni_ray = sign(bidray.B[index_non_zero, k+n]) * bidray.B[index_non_zero, :]

            # Update bidirectional rays
            new_bidray_mat = remove_row(bidray.B, index_non_zero)
            for i in 1:bidray.rows-1
                if new_bidray_mat[i,n+k] != 0
                    new_bidray_mat[i, :] = combine(k+n, vec(new_uni_ray), vec(new_bidray_mat[i,:]))
                end
            end
            bidray.B = new_bidray_mat
            bidray.rows -= 1

            # Update unidirectional rays
            new_uniray_mat = add_row(uniray.B, new_uni_ray)
            for i in 1:uniray.rows
                if uniray.B[i, n+k] != 0
                    new_uniray_mat[i,:] = combine(new_uniray_mat, k+n, uniray.rows+1, i)
                end
            end
            uniray.B = new_uniray_mat
            uniray.rows += 1
        else
            if verbosity > 0
                println("All bidirectional rays are saturated by constraint $(k)")
            end
            
            pos, nil, neg = get_pos_neg_indices(uniray, k+n)
            if verbosity > 1
                println("pos, nil, neg: $(pos), $(nil), $(neg)")
            end

            if length(neg) == 0
                if verbosity > 0
                    println("Constraint $k is violated by no bidirectional or unidirectional rays so is redundant")
                end
                push!(redundant_rows, k)
                continue
            end

            if length(pos) + length(nil) == 0
                if verbosity > 0
                    println("Constraint $k violates all bidirectional and unidirectional rays so Cone is zero set")
                end
                return zeros(n)
            end
            
            if verbosity > 0
                println("Adding constraint $(k) which is currently violated by $(length(neg)) of $(uniray.rows) extremal rays")
                println("New unidirectional basis will consist of a maximum of $(uniray.rows) - $(length(neg)) + $(length(pos)) * $(length(neg)) = $(uniray.rows - length(neg) + length(pos)*length(neg)) extremal rays")
            end;
            
            max_rows = length(pos) + length(nil) + length(pos) * length(neg)
            new_B = Array(T, max_rows, uniray.m + uniray.n)


            # Handle case of two rows separately
            if uniray.rows == 2
                if verbosity > 0; println("Handling two row case..."); end;
                new_B[1,:] = uniray.B[pos[1],:]
                new_row = combine(uniray.B, k+n, pos[1], neg[1])
                if verbosity > 1; println("Candidate row: ", new_row); end;
                if !collinear(vec(uniray.B[pos[1],:]), vec(new_row))
                    new_B[2,:] = new_row
                    uniray.rows = 2
                    uniray.B = new_B
                else
                    uniray.rows = 1
                    uniray.B = new_B[1,:]
                end
                continue
            end
            
            # Add rows with non-negative intersections
            # with leading column to next matrix
            for (l,i) in enumerate([pos,nil]); new_B[l,:] = uniray.B[i,:]; end;
            new_rows = length(pos) + length(nil)

            # Go through pairs of row with positive and negative
            # coefficients in leading column and combine if necessary
            for i1 in pos, i2 in neg
                if combinable(uniray, i1, i2, verbosity - 1)
                    new_rows += 1
                    new_B[new_rows, :] = combine(uniray.B, k+n, i1, i2)
                end
            end
            uniray.B = new_B[1:new_rows,:]
            uniray.rows = new_rows

            if verbosity > 0
                @printf("Added %d / %d positive combinations of extremal rays\n",
                        new_rows - length(pos) - length(nil), length(pos) * length(neg))
                println("------------------------------------------------")
            end
        end
        union!(redundant_rows, redundant_constraint_check(uniray,k, verbosity))
    end
    if verbosity > 0
        println("Result")
        println("---------------")
        println("Bidirectional rays")
        pprint(bidray)
        println("Unidirectional rays")
        pprint(uniray)
    end
    return norm_cols(LHS(bidray)'), norm_cols(LHS(uniray)'), redundant_rows
end
