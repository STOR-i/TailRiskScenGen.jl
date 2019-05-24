"""Reduces a scenario set by aggregating all scenarios outside a specified risk region"""
function aggregate_scenarios(scenarios::Matrix{Float64}, Ω::AbstractRiskRegion)
    num_risk = 0
    num_non_risk = 0
    dim, num_scen = size(scenarios)
    new_scenarios = Array{Float64}(undef, dim, num_scen)
    non_risk_sum = zeros(dim)
    for s in 1:num_scen
        if scenarios[:,s] ∈ Ω
            new_scenarios[:,(num_risk+=1)] = scenarios[:,s]
        else
            non_risk_sum += scenarios[:,s]
            num_non_risk += 1
        end
    end
    if num_non_risk == 0
        warn("No scenarios in non-risk region to aggregate")
        new_probs = fill(1.0/(num_scen), num_scen)
        return new_scenarios, new_probs
    else
        new_num_scen = num_risk + 1
        new_scenarios[:,new_num_scen] = non_risk_sum/num_non_risk
        new_probs = fill(1.0/(num_risk + num_non_risk), new_num_scen)
        new_probs[new_num_scen] = num_non_risk/(num_risk + num_non_risk)
        # Resize array to have appropriate number of scenarios - pointer hack causes bugs - just copy into new array
        # new_scenarios = pointer_to_array(pointer(new_scenarios), (dim, new_num_scen))
        new_scenarios = new_scenarios[:, 1:new_num_scen]
        return new_scenarios, new_probs
    end
end

function _aggregation_sampling(dist::Sampleable{Multivariate, Continuous}, Ω::AbstractRiskRegion, num_scen::Int64)
    dim = length(dist)
    scenarios = Array{Float64}(undef, dim, num_scen)
    non_risk_sum = fill(0.0, dim)
    num_risk = 0
    num_non_risk = 0
    while num_risk + 1 < num_scen
        scenarios[:, num_risk + 1] = rand(dist)
        if scenarios[:, num_risk+1] ∈ Ω
            num_risk += 1
        else
            non_risk_sum += scenarios[:, num_risk + 1]
            num_non_risk += 1
        end
    end
    scenarios[:,num_scen] = non_risk_sum/num_non_risk
    probs = fill(1.0/(num_risk + num_non_risk), num_scen)
    probs[num_scen] = num_non_risk/(num_risk + num_non_risk)
    results = Dict()
    results["scenarios"] = scenarios
    results["probs"] = probs
    results["num_sampled"] = num_risk + num_non_risk
    results
end

"""Constructs scenarios via aggregation sampling for a given distribution and risk region"""
function aggregation_sampling(dist::Sampleable{Multivariate, Continuous}, Ω::AbstractRiskRegion, num_scen::Int64)
    results = _aggregation_sampling(dist, Ω, num_scen)
    return results["scenarios"], results["probs"]
end

# function _parallel_aggregation_sampling(dist::Sampleable{Multivariate, Continuous},
#                                         Ω::AbstractRiskRegion, num_scen::Int,
#                                         num_threads::Int)
#     d = length(dist)
#     n, r = div(num_scen, num_threads), num_scen % num_threads
#     threads = []
#     push!(threads, @spawnat 1 _aggregation_sampling(dist, Ω, n+r))
#     for i in 1:num_threads
#         push!(threads, @spawnat i _aggregation_sampling(dist, Ω, n+1))
#     end

#     results_list = [fetch(p) for p in threads]
#     total_sampled = sum(res["num_sampled"] for res in results_list)
#     scenarios = Array{Float64}(undef, d, num_scen)
#     scenarios[:, 1:(n+r-1)] = results_list[1]["scenarios"][:, 1:(n+r-1)]
#     total_non_risk = results_list[1]["num_sampled"] - (n+r-1)
#     scenarios[:, num_scen] = total_non_risk*results_list[1]["scenarios"][:, n+r]
#     c = n+r
#     for i in 2:num_threads
#         scenarios[:,c:(c+n-1)] = results_list[i]["scenarios"][:,1:n]
#         num_non_risk = results_list[i]["num_sampled"] - (n+r-1)
#         total_non_risk += num_non_risk
#         scenarios[:, num_scen] += num_non_risk*results_list[i]["scenarios"][:, n+1]
#     end

#     scenarios[:, num_scen]/=total_sampled
#     probs = fill(1.0/total_sampled, num_scen)
#     probs[num_scen] = total_non_risk/total_sampled
#     results = Dict()
#     results["num_sampled"] = total_sampled
#     results["scenarios"] = scenarios
#     results["probs"] = probs
#     return results
# end

