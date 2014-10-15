using Cpp
using Distributions
using Gadfly

function aggregate_scenarios(scenarios::Array{Float64, 2}, probs::Array{Float64, 1},
                             mean::Array{Float64, 1}, cov::Array{Float64, 2},
                             cone_A::Array{Float64, 2}, alpha::Float64)
    dim = size(mean,1)
    num_scen = size(scenarios, 2)
    num_gen = size(cone_A, 2)
    new_scenarios = Array(Float64, dim, num_scen)
    new_probs = Array(Float64, num_scen)
    new_num_scen = @cpp ccall( ("aggregation_julia", "bin/libScenGen.so"),
                              Int,
                              (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                               Int, Int, Int,
                               Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Float64),
                              new_scenarios, new_probs, scenarios, probs, dim, num_scen, num_gen,
                              mean, cov, cone_A, alpha)
    return new_scenarios[:, 1:new_num_scen], new_probs[1:new_num_scen]
end

mu = [0.0, 0.5]
cov = [[1 0.5], [0.5 1]]
num_scen = 1000

normal = MvNormal(mu, cov)
scenarios = rand(normal, num_scen)
probs = fill(1.0/num_scen, num_scen)

new_scen, new_probs = aggregate_scenarios(scenarios, probs, mu, cov, eye(2), 1.96)
x_min, y_min = minimum(scenarios,2)
x_max, y_max = maximum(scenarios,2)
x_scale = Scale.x_continuous(minvalue=x_min, maxvalue=x_max)
y_scale = Scale.y_continuous(minvalue=y_min, maxvalue=y_max)

p1 = plot(x=scenarios[1,:], y=scenarios[2,:], x_scale, y_scale, Guide.title("Full Scenario Set"))
p2 = plot(x=new_scen[1,:], y=new_scen[2,:], x_scale, y_scale, Guide.title("Aggregated Scenario Set"))
draw(PDF("data.pdf", 6inch, 6inch), vstack(p1,p2))

