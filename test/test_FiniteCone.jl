# Tests whether points project correctly

print("\tTest 1 - #generators == dim == 2...")
cone_1 = FiniteCone([[0.0, 1.0] [1.0, 1.0]])
test_data_1 = [([-1.0, 3.0], [0.0, 3.0]),
               ([-1.0, -1.0], [0.0, 0.0]),
               ([2.0, 0.0], [1.0, 1.0]),
               ([0.25, 1.0], [0.25, 1.0])]
for (x, proj_x) in test_data_1
    @test_approx_eq project(cone_1, x) proj_x
end
print("success!\n")

print("\tTest 2 - 3 = #generators > dim = 2...")
cone_2 = FiniteCone([[0.0, 1.0] [1.0, 1.0] [0.5, 1.0]])
test_data_2 = [([-1.0, 3.0], [0.0, 3.0]),
               ([-1.0, -1.0], [0.0, 0.0]),
               ([2.0, 0.0], [1.0, 1.0]),
               ([0.25, 1.0], [0.25, 1.0])]
for (x, proj_x) in test_data_2
    @test_approx_eq project(cone_2, x) proj_x
end
print("success!\n")
    

print("\tTest 3 - 2 = #generators < dim = 3...")
cone_3 = FiniteCone([[0.0, 0.0, 1.0] [1.0, 1.0, 1.0]])
test_data_3 = [([-1.0, -1.0, 1.0], [0.0, 0.0, 1.0]),
               ([3.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
               ([-1.0, -1.0, -1.0], [0.0, 0.0, 0.0]),
               ([1.0, 2.0, 3.0], [1.5, 1.5, 3.0]),
               ([2.5, 2.5, 6.0], [2.5, 2.5, 6.0])]
for (x, proj_x) in test_data_3
    @test_approx_eq project(cone_3, x) proj_x
end
print("success!\n")
