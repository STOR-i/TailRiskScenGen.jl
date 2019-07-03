# TailRiskScenGen.jl

This Julia package provides functionality for generating scenarios for stochastic programs which use a tail risk measure such as the conditional value-at-risk. For problems involving tail risk measure, the value of the tail risk measure depends only a subset of the support of the distribution called the **risk region**. By prioritising the generation of scenarios in the risk region, one can approximate much better the value of a tail risk measure in a stochastic program.

This package allows users to construct risk regions and use these scenario generaion and reduction.
Currently risk regions are provided for two types of problems:

- Portfolio selection problem with Elliptical return distributions
- Stochastic programs with monotonic loss/recourse functions.

In the former case it has been shown that the methodology still performs well when the returns have near-elliptical distributions.

For more details on this methodology see the following two papers:

- [Problem-driven scenario generation: an analytical approach for stochastic programs with tail risk measure](http://arxiv.org/abs/1511.03074)
- [Scenario Generation for Single-Period Portfolio Selection Problems with Tail Risk Measures: Coping with High Dimensions and Integer Variables](https://pubsonline.informs.org/doi/10.1287/ijoc.2017.0790)

# Risk regions

The risk region for our class of problems have the following form:

![Basic form of risk regions](/docs/riskregion1.png?raw=true "Formula for risk regions")

Given this formulation It is difficult to test whether a scenario `y` is in the risk region. When `K` is a convex cone, the above expression can be rewritten as follows:

![Nice form of risk regions](/docs/riskregion2.png?raw=true "Nice formula for risk regions")

where `K' = PK` and `p_k(y)` is the projection of the point y onto the cone K.

The construction of a risk region requires four parameters: a vector μ, a positive definite matrix Σ=PᵀP, a convex cone K and positive scalar α. The package currently contains two types of cone objects: a FiniteCone and a PolyhedralCone. A FiniteCone object is constructed from a matrix of cone generators where each column corresponds to a cone generator.

```
K = FiniteCone(eye(2)) # Constructs cone of positive quadrant in R^2
```

We now construct a RiskRegion for Normally distributed returns at the level β = 0.95:

```
μ = [0.05, 0.02]
Σ = [[0.5, 0.2] [0.2, 0.5]]
β = 0.95
Ω = RiskRegion(μ, Σ, K, quantile(Normal(), β))
```

Convenience constructors are provided for multivariate normal and t distributions:

```
Y = MvNormal(μ, Σ)
Ω = RiskRegion(Y, K, β)
```
# Scenario generation and scenario reduction

Scenario reduction requires a matrix of (equiprobable) scenarios, where each column
corresponds to one scenario, and a risk region. The aggregate_scenarios function
outputs a scenario matrix where all scenarios in the non-risk region have been aggregated into
a single point, and the corresponding vector of scenario probabilities.

```
n = 1000                                # Initial number of scenarios
scenarios = rand(Y, num_s)                         # Generate a matrix of scenarios
reduced_set, reduced_probs = aggregate_scenarios(scenarios, Ω) # Reduce scenario set
```

A scenario set of specified size can be generated via the **aggregation sampling** algorithm.
This takes a risk region and a distribution, and samples points from this, stores any point in the risk region
and aggregates any in the non-risk region. The algorithm terminates when the required number of points in the risk region have been sampled.

```
n = 500 # Scenario set size
scenarios, probs = aggregation_sampling(Y, Ω, n)
```