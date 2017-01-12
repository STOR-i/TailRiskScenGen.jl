function vec_isless(x::AbstractArray, y::AbstractArray)
    for (xi,yi) in zip(x,y)
        if xi >= yi return false end
    end
    return true
end

type SurvivorApproximator
    sample::MatF64
    α::Float64
end

function (surv::SurvivorApproximator)(x::VecF64)
    count=0
    N = size(surv.sample, 2)
    # count = @parallel (+) for i in 1:N
    #     convert(Int, vec_isless(x, view(surv.sample, :, i)))
    # end
    for i in 1:N
        count += convert(Int, vec_isless(x, view(surv.sample, :, i)))
    end
    prob=count/N
    se=sqrt((prob-prob^2)/N)
    err=norminvcdf((1 + surv.α)/2)*se
    return prob,err
end

function below_nonrisk_frontier{T<:Real}(x::AbstractVector{T}, frontier::LinkedList{VecF64})
    for p in frontier
        if vec_isless(x,p) return true end
    end
    return false
end

function above_risk_frontier{T<:Real}(x::AbstractVector{T}, frontier::LinkedList{VecF64})
    for p in frontier
        if vec_isless(p, x) return true end
    end
    return false
end


function decapitate!(l::Cons)
    l.head=l.tail.head
    l.tail=l.tail.tail
    return l
end
