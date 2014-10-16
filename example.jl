using EllipticalScenGen
using Distributions
## using Gadfly


mu = [0.0, 0.5]
cov = [[1 0.5], [0.5 1]]
num_scen = 1000

normal = MvNormal(mu, cov)
scenarios = rand(normal, num_scen)
probs = fill(1.0/num_scen, num_scen)

# print(scenarios)

new_scen, new_probs = aggregate_scenarios(scenarios, probs, mu, cov, eye(2), quantile(Normal(), 0.95))
# print(new_scenarios)


## x_min, y_min = minimum(scenarios,2)
## x_max, y_max = maximum(scenarios,2)
## x_scale = Scale.x_continuous(minvalue=x_min, maxvalue=x_max)
## y_scale = Scale.y_continuous(minvalue=y_min, maxvalue=y_max)
## p1 = plot(x=scenarios[1,:], y=scenarios[2,:], x_scale, y_scale, Guide.title("Full Scenario Set"))
## p2 = plot(x=new_scen[1,:], y=new_scen[2,:], x_scale, y_scale, Guide.title("Aggregated Scenario Set"))
## draw(PDF("data.pdf", 6inch, 6inch), vstack(p1,p2))


