using Cpp
using Distributions

lib_path =  joinpath(dirname(@__FILE__()), "../bin/libScenGen.so")
const lib = normpath(lib_path)

function aggregate_scenarios(scenarios::Array{Float64, 2}, probs::Array{Float64, 1},
                             mean::Array{Float64, 1}, cov::Array{Float64, 2},
                             cone_A::Array{Float64, 2}, alpha::Float64)
    dim = size(mean,1)
    num_scen = size(scenarios, 2)
    num_gen = size(cone_A, 2)
    new_scenarios = Array(Float64, dim, num_scen)
    new_probs = Array(Float64, num_scen)
    new_num_scen = @cpp ccall( ("aggregation_julia", lib),
                              Int,
                              (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                               Int, Int, Int,
                               Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Float64),
                              new_scenarios, new_probs, scenarios, probs, dim, num_scen, num_gen,
                              mean, cov, cone_A, alpha)
    return new_scenarios[:, 1:new_num_scen], new_probs[1:new_num_scen]
end

function aggregate_scenarios(scenarios::Matrix{Float64}, Ω::RiskRegion)
    num_risk::Int64 = 0
    num_non_risk::Int64 = 0
    dim, num_scen = size(scenarios)
    new_scenarios = Array(Float64, dim, num_scen)
    non_risk_sum = fill(0.0, dim)
    for s in 1:num_scen
        if in_RiskRegion(Ω, scenarios[:,s])
            new_scenarios[:,(num_risk+=1)] = scenarios[:,s]
        else
            non_risk_sum += scenarios[:,s]
            num_non_risk += 1
        end
    end
    new_num_scen = num_risk + 1
    new_scenarios[:,new_num_scen] = non_risk_sum/num_non_risk
    new_probs = fill(1.0/(num_risk + num_non_risk), new_num_scen)
    new_probs[new_num_scen] = num_non_risk/(num_risk + num_non_risk)
    # Resize array to have appropriate number of scenarios - uses pointer hack
    new_scenarios = pointer_to_array(pointer(new_scenarios), (dim, new_num_scen))
    new_scenarios, new_probs
end

function aggregation_sampling(dist::Sampleable{Multivariate, Continuous}, Ω::RiskRegion, num_scen::Int64)
    dim = length(dist)
    scenarios = Array(Float64, dim, num_scen)

    non_risk_sum = Array(Float64, dim)
    num_risk = 0
    num_non_risk = 0
    while num_risk + 1 < num_scen
        scenarios[:, num_risk + 1] = rand(dist)
        if in_RiskRegion(Ω, scenarios[:, num_risk + 1])
            num_risk += 1
        else
            non_risk_sum += scenarios[:, num_risk + 1]
            num_non_risk += 1
        end
    end
    scenarios[:,num_scen] = non_risk_sum/num_non_risk
    probs = fill(1.0/(num_risk + num_non_risk), num_scen)
    probs[num_scen] = num_non_risk/(num_risk + num_non_risk)
    scenarios, probs
end
