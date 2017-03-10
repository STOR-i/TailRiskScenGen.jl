type MonotonicRiskRegion <: AbstractRiskRegion
    β::Float64
    nonrisk_frontier::LinkedList{VecF64}
    risk_frontier::LinkedList{VecF64}
    surv::SurvivorApproximator
    function MonotonicRiskRegion(surv::SurvivorApproximator, β::Float64)
        return new(β,nil(VecF64),nil(VecF64),surv)
    end
end

function add_to_risk_frontier!{T<:Real}(Ω::AbstractRiskRegion, x::AbstractVector{T})
    l = Ω.risk_frontier
    if isa(l, Nil{VecF64})
        return Ω.risk_frontier=list(x[:])
    end

    while !isa(l, Nil{VecF64})
        pf = head(l)
        if vec_isless(x, pf)
            if isa(tail(l),Nil{VecF64})
                return Ω.risk_frontier=list(x[:])
            else
                decapitate!(l)
            end
        end
        l = tail(l)
    end
    return Ω.risk_frontier = cons(x[:], Ω.risk_frontier)
end

function add_to_nonrisk_frontier!{T<:Real}(Ω::AbstractRiskRegion, x::AbstractVector{T})
    l = Ω.nonrisk_frontier
    if isa(l, Nil{VecF64})
        return Ω.nonrisk_frontier=list(x[:])
    end
    
    while !isa(l, Nil{VecF64})
        pf = head(l)
        if vec_isless(pf, x)
            if isa(tail(l),Nil{VecF64})
                return Ω.nonrisk_frontier=list(x[:])
            else
                decapitate!(l)
            end
        end
        l = tail(l)
    end
    return Ω.nonrisk_frontier = cons(x[:], Ω.nonrisk_frontier)
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

function prob_nonrisk(dist::Distribution{Multivariate, Continuous}, β_list::Vector{Float64}, num_points_surv::Int = 50000, num_points_nonrisk::Int = 10000, α::Float64=0.99)
    sample = rand(dist, num_points_nonrisk)
    surv = SurvivorApproximator(rand(dist, num_points_surv), α)
    Ω = Dict{Float64, MonotonicRiskRegion}()
    probs = Dict{Float64,Float64}()
    
    for β in β_list
        Ω[β] = MonotonicRiskRegion(surv, β)
        probs[β] = 0.0
    end

    for i in 1:num_points_nonrisk
        point = view(sample, :, i)
        calc_surv_prob = false
        surv_prob = 0.0
        for β in β_list
            if calc_surv_prob
                probs[β] += (surv_prob > 1-β ? 1 : 0)
                continue
            else
                shortcut = false
                # Does sampled point dominate any point in risk region?
                if below_nonrisk_frontier(point, Ω[β].nonrisk_frontier)
                    shortcut=true
                    probs[β] += 1
                elseif above_risk_frontier(point, Ω[β].risk_frontier)
                    shortcut=true
                end

                if !shortcut
                    surv_prob, err = surv(point)
                    if surv_prob > 1 - β
                        probs[β]+=1
                        if surv_prob - err > 1-β
                            add_to_nonrisk_frontier!(Ω[β], point)
                        end
                    else
                        if surv_prob + err < 1-β
                            add_to_risk_frontier!(Ω[β], point)
                        end
                    end
                    calc_surv_prob=true
                end
            end
        end
    end
    
    for β in β_list
        probs[β]/= num_points_nonrisk
    end
    return probs
end
