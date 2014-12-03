lib_path =  joinpath(dirname(@__FILE__()), "../bin/libScenGen.so")
const lib = normpath(lib_path)

function lcp_solve(M::Matrix{Float64}, q::Vector{Float64})
    n = length(q)
    w,z = Array(Float64,n), Array(Float64, n)
    success = ccall((:lcp_julia_solve, lib), Int, (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int), w, z, M, q, n)
    if success != 1
        print("WARNING: lcp_solve did not succeed\n")
    end
    return w,z
end
