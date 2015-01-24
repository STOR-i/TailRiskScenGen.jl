using Distributions

A = [[3.0 -4.0 1.0 0.0],
     [2.0 0.0 4.0 -1.0],
     [-4.0 -7.0 2.0 4.0],
     [-1.0 0.0 20.0 2.0],
     [6.0 -5.0 -4.0 2.0]]

X = chernikova(A)

# Function which checks anything generated from Chernikova output is in Polyhedral cone
function test_polyhedral_contains_finite_cone(A::Matrix{Float64}, num_points::Int)
    X = chernikova(A)
    gens, dim = size(X)
    dist = Uniform(0.0, 100.0)
    for i in 1:num_points
        λ = rand(dist, gens)
        @test all(A*X'λ .>= 0.0)
    end
end


test_polyhedral_contains_finite_cone(A, 10000)
