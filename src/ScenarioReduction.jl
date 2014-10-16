using Cpp
using Distributions

dir = dirname(@__FILE__())
lib =  joinpath(dir, "../bin/libScenGen.so")
lib = normpath(lib)

function aggregate_scenarios(scenarios::Array{Float64, 2}, probs::Array{Float64, 1},
                             mean::Array{Float64, 1}, cov::Array{Float64, 2},
                             cone_A::Array{Float64, 2}, alpha::Float64)
    dim = size(mean,1)
    num_scen = size(scenarios, 2)
    num_gen = size(cone_A, 2)
    new_scenarios = Array(Float64, dim, num_scen)
    new_probs = Array(Float64, num_scen)
    new_num_scen = @eval @cpp ccall( ("aggregation_julia", $lib),
                                    Int,
                                    (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                                     Int, Int, Int,
                                     Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Float64),
                                    $new_scenarios, $new_probs, $scenarios, $probs, $dim, $num_scen, $num_gen,
                                    $mean, $cov, $cone_A, $alpha)
    return new_scenarios[:, 1:new_num_scen], new_probs[1:new_num_scen]
end
