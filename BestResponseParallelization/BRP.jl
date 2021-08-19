using OneTest
using Random
using DataFrames, CSV

# This brief experiment suggests that the best-response
# calculations parallelize well.

function timetester(n=500)
    N = Int64[]
    T_par = Float64[]

    for t in 1:n
        C = rand(20:150)
        m = Market(C)
        s = randexp(C)
        p = rand(C)

        push!(N, C)
        push!(T_par, minimum(@elapsed bestresponse_it(m, s, p) for i in 1:4))
    end

    return N, T_par
end

nt = Threads.nthreads()
@show nt
N, T_par = timetester()

CSV.write("./BestResponseParallelization/times/$(nt)threads.csv", DataFrame([:C => N, :time => T_par]))
