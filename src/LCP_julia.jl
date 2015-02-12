lib_path =  joinpath(dirname(@__FILE__()), "../bin/libScenGen.so")
const lib = normpath(lib_path)

function lcp_solve(M::Matrix{Float64}, q::Vector{Float64})
    n = length(q)
    w,z = Array(Float64,n), Array(Float64, n)
    return_code = ccall((:lcp_julia_solve, lib), Int, (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int), w, z, M, q, n)
    if return_code == 1
        warn("lcp_solve reached maximum number of iterations")
    elseif return_code == 2
        warn("lcp_solve failed due to negative diagonal term")
    end
    return w,z
end
