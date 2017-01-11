type MonotonicRiskRegion
    β::Float64
    nonrisk_frontier::LinkedList{VecF64}
    risk_frontier::LinkedList{VecF64}
    surv::SurvivorApproximator
    function MonotonicRiskRegion(surv::SurvivorApproximator, β::Float64)
        return new(β,nil(VecF64),nil(VecF64),surv)
    end
end

function add_to_risk_frontier!(Ω::MonotonicRiskRegion, x::VecF64)
    l = Ω.risk_frontier
    if isa(l, Nil{VecF64})
        return Ω.risk_frontier=list(x)
    end

    while !isa(l, Nil{VecF64})
        pf = head(l)
        if vec_isless(x, pf)
            if isa(tail(l),Nil{VecF64})
                return Ω.risk_frontier=list(x)
            else
                decapitate!(l)
            end
        end
        l = tail(l)
    end
    return Ω.risk_frontier = cons(x, Ω.risk_frontier)
end

function add_to_nonrisk_frontier!(Ω::MonotonicRiskRegion, x::VecF64)
    l = Ω.nonrisk_frontier
    if isa(l, Nil{VecF64})
        return Ω.nonrisk_frontier=list(x)
    end
    
    while !isa(l, Nil{VecF64})
        pf = head(l)
        if vec_isless(pf, x)
            if isa(tail(l),Nil{VecF64})
                return Ω.nonrisk_frontier=list(x)
            else
                decapitate!(l)
            end
        end
        l = tail(l)
    end
    return Ω.nonrisk_frontier = cons(x, Ω.nonrisk_frontier)
end

function ∈(x::VecF64, Ω::MonotonicRiskRegion)
    if below_nonrisk_frontier(x, Ω.nonrisk_frontier) return false
    elseif above_risk_frontier(x, Ω.risk_frontier) return true
    else
        prob, err = Ω.surv(x)
        if prob > 1 - Ω.β
            if prob - err > (1 - Ω.β)
                add_to_nonrisk_frontier!(Ω, x)
            end
            return false
        else
            if prob + err < (1-Ω.β)
                add_to_risk_frontier!(Ω, x)
            end
            return true
        end
    end
end
