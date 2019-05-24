function nonrisk_clustering(dist::Sampleable{Multivariate, Continuous}, Ω::AbstractRiskRegion,
                            num_risk::Int64, num_non_risk::Int64)
    max_non_risk = (num_risk + num_non_risk) * 10
    dim = length(dist)
    scenarios = Array{Float64}(undef, dim, num_risk + num_non_risk)
    non_risk_scen = Array{Float64}(undef, dim, max_non_risk)
    r = 0                  # Counter for risk scenarios
    nr = 0   # Counter for non-risk scenarios

    while r < num_risk
        scen = rand(dist)
        if scen ∈ Ω
            scenarios[:,r+1] = scen
            r += 1
        else
            non_risk_scen[:, nr+1] = scen
            nr += 1
        end
    end

    if nr <= num_non_risk
        scenarios[:,r+1:r+nr] = non_risk_scen[:,1:nr]
        probs = fill(1.0/(r+nr), r+nr)
        return scenarios, probs
    else
        num_clusters = num_non_risk
        cluster_results = kmeans(non_risk_scen[:,1:nr], num_clusters)
        scenarios[:, (num_risk+1):(num_risk+num_non_risk)] = cluster_results.centers
        probs = Array{Float64}(undef, num_risk + num_clusters)
        probs[1:r] .= 1.0/(r + nr)
        probs[r+1:r + num_clusters] .= cluster_results.cweights/(r + nr)
    end
    return scenarios, probs
end

function nonrisk_clustering(scenarios::Matrix{Float64}, Ω::AbstractRiskRegion, num_non_risk::Int64)
    dim, num_scen = size(scenarios)
    new_scenarios = Array{Float64}(undef, dim, num_scen)
    non_risk_scenarios = Array{Float64}(undef, dim, num_scen)
    r = 0   # Risk scenario counter
    nr = 0  # Non-risk scenario counter
    for s in 1:num_scen
        scen = scenarios[:,s]
        if scen ∈ Ω
            new_scenarios[:,r+1] = scen
            r += 1
        else
            non_risk_scenarios[:, nr+1] = scen
            nr += 1
        end
    end
    non_risk_scenarios = non_risk_scenarios[:, 1:nr]
    if nr > num_non_risk
        cluster_results = kmeans(non_risk_scenarios, num_non_risk)
        new_scenarios[:,(r+1):(r+num_non_risk)] = cluster_results.centers
        # The following line causes a garbage-collection related bug, so we use something less efficient...
        # new_scenarios = pointer_to_array(pointer(new_scenarios), (dim, r + num_non_risk))
        new_scenarios = new_scenarios[:, 1:r + num_non_risk]
        new_probs = Array{Float64}(undef, r + num_non_risk)
        new_probs[1:r] = fill(1.0/(r + nr), r)
        new_probs[r+1:r + num_non_risk] = cluster_results.cweights/(r + nr)
    else
        new_scenarios[:,(r+1):(r+nr)] = non_risk_scenarios
        # new_scenarios = pointer_to_array(pointer(new_scenarios), (dim, r + nr)) # Bug
        new_scenarios = new_scenarios[:, 1:(r+nr)]
        new_probs = fill(1.0/(r+nr))
    end
    return new_scenarios, new_probs
end
