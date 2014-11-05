using Clustering

function nonrisk_clustering(dist::Sampleable{Multivariate, Continuous}, Ω::RiskRegion,
                            num_risk::Int64, num_non_risk::Int64)
    max_non_risk = (num_risk + num_non_risk) * 10
    dim = length(dist)
    scenarios = Array(Float64, dim, num_risk + num_non_risk)
    non_risk_scen = Array(Float64, dim, max_non_risk)
    r = 0                  # Counter for risk scenarios
    nr = 0   # Counter for non-risk scenarios

    while r < num_risk
        scen = rand(dist)
        if in_RiskRegion(Ω, scen)
            scenarios[:,r+1] = scen
            r += 1
        else
            non_risk_scen[:, nr+1] = scen
            nr += 1
        end
    end

    num_clusters = min(num_non_risk, nr)
    cluster_results = kmeans(non_risk_scen[:,1:nr], num_clusters)
    scenarios[:, (num_risk+1):(num_risk+num_non_risk)] = cluster_results.centers
    probs = Array(Float64, num_risk + num_clusters)
    probs[1:num_risk] = 1.0/(r + nr)
    probs[num_risk+1:num_risk + num_clusters] = cluster_results.cweights/(num_risk + nr)
    return scenarios, probs
end
