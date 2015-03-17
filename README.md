# EllipticalScenGen

This Julia package provides tools for generating scenario sets. The methodology employed is in particular adapted to portfolio selection problems using tail risk measures where the returns have near-elliptical distributions. For this problems involving tail risk measure, the value of the tail risk measure depends only a subset of the support of the distribution called the **risk region**. This package exploits this property by providing functions which prioritise the generation of scenarios in these risk regions for the class of problems described above.

The author can provide a draft paper with more details.

# Risk regions

The risk region for our class of problems have the following form:

![Basic form of risk regions](/docs/riskregion1.pdf?raw=true "Formula for risk regions")

Given this formulation It is difficult to test whether a scenario `y` is in the risk region. When `K` is a convex cone, the above expression can be rewritten as follows:

![Nice form of risk regions](/docs/riskregion2.pdf?raw=true "Nice formula for risk regions")

# Acknowledgements

This project makes uses code from the Siconos Numerics project
(http://siconos.gforge.inria.fr/) developed by INRIA for solving
linear complementarity complementarity problems using Lemke's
algorithm.