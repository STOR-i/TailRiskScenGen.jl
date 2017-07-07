using TailRiskScenGen
using Distributions

# Inputs - change as necessary
dim = 20
gens = 30
samples = 100000
srand(1)

# Randomly generate distribution and cone
Z = Normal()
μ = rand(Z, dim)
A = rand(Z, (dim,dim))
Σ = A'A
dist = MvNormal(μ, Σ)
K = FiniteCone(rand(Z, (dim, gens)))

function test(dist::MvNormal, K::FiniteCone, samples::Int64)
    dim = K.dim
    for i in 1:samples
        project(K, rand(dim))
    end
end    

# Small test to compile necessary functions
test(dist, K, 5)

# Real test
@time test(dist, K, samples)




