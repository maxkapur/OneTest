module OneTest

import Base: length
import Base: sort
using Random, Distributions
using LinearAlgebra

export Market, demand, appeal, demandmatrix, tatonnement, equilibrium, makepreflists
export utility, incentivegradient, localequilibriumsearch, bestresponse_it, bestresponse, sigmainvopt


"""
Contains static information about a school-choice market.
"""
struct Market{T<:AbstractFloat}
    qualities ::Array{T, 1}    # δ
    capacities::Array{T, 1}    # q
    gamma     ::Array{T, 1}    # γ
end


length(market::Market) = length(market.qualities)


"""
Sort the market by gamma over capacities, i.e. in 
ascending order of selectivity at equilibrium.
"""
function sort(market::Market)
    sort_order = sortperm(market.gamma ./ market.capacities)
    return Market(market.qualities[sort_order],
                  market.capacities[sort_order],
                  market.gamma[sort_order])
end


"""
    Market(qualities, capacities)

Instantiate a school-choice market. `gamma` is normalized to sum to one.
"""
function Market(qualities::Array{T, 1}, capacities::Array{T, 1}) where {T<:AbstractFloat}
    length(qualities) == length(capacities) || return throw(DimensionMismatch)
    gamma = exp.(qualities)
    return Market{T}(qualities, capacities, gamma)
end
    
    
"""
    Market(m)

Generate a random market with m schools.
"""
function Market(m::Int) 
    qualities = rand(m)
    capacities = 2 * rand(m)/m
    return Market(qualities, capacities)
end


"""
    demand(market, cutoffs)

Return demand for each school given a set of cutoffs and ignoring capacity, using
multinomial logit choice model with one student profile and a single test score.
"""
function demand(market      ::Market,
                cutoffs     ::AbstractArray{<:AbstractFloat, 1};
                )::AbstractArray{<:AbstractFloat, 1}
    
    m = length(market)
    demands = zeros(m)

    sort_order = sortperm(cutoffs)
    cutoffs[sort_order]

    γ = exp.(market.qualities)
    demands = zeros(m)

    consideration_set_probabilities = diff([cutoffs[sort_order]; 1])

    for c in 1:m, d in c:m     # For each score threshold
        demands[sort_order[c]] += consideration_set_probabilities[d] *
                                  γ[sort_order[c]] / sum(γ[sort_order[1:d]])
    end

    return demands
end


"""
    appeal(market, cutoffs, sigma=1)

Return appeal of entering class at each school given a set of cutoffs, using
multinomial logit choice model with one student profile and a single test score.

Sigma encodes the school's valuation function over the interval `[0, 1]`, namely
``v(x) = 1 + x^ σ``. 
"""
function appeal(market      ::Market,
                cutoffs     ::AbstractArray{<:AbstractFloat, 1};
                sigma=1     ::Union{<:Integer, <:AbstractFloat, AbstractArray{<:AbstractFloat, 1}}
                )::AbstractArray{<:AbstractFloat, 1}
    
    m = length(market)
    demands = zeros(m)

    sort_order = sortperm(cutoffs)
    cutoffs[sort_order]

    γ = exp.(market.qualities)
    appeals = zeros(m)

    diff_of_squares = diff([cutoffs[sort_order]; 1] .^ (sigma .+ 1))

    for c in 1:m, d in c:m     # For each score threshold
        appeals[sort_order[c]] += diff_of_squares[d] *
                                  γ[sort_order[c]] / sum(γ[sort_order[1:d]])
    end

    return appeals ./ (sigma .+ 1)
end


"""
    demandmatrix(market, cutoffs)

Returns the matrices `A` and the permutation `sort_order` such that

````
demand(market, p)[sort_order] == A * cutoffs[sort_order] + market.gamma[sort_order]/sum(market.gamma)
````

and

````
appeal(market, p)[sort_order] == (A * cutoffs[sort_order] .^2 + market.gamma[sort_order]/sum(market.gamma))/2
````

The argument `cutoffs` is used only to determine the sort order, so it may be
replaced with any vector whose entries sort the same way.
"""
function demandmatrix(market, cutoffs)
    m = length(market)
    
    sort_order = sortperm(cutoffs)
    
    A = UpperTriangular(zeros(m, m))
    
    for c in 1:m
        A[c, c] = -market.gamma[sort_order[c]] / sum(market.gamma[sort_order[1:c]])
    end
    
    for i in 1:m, j in i+1:m
        A[i, j] = market.gamma[sort_order[i]] *
                    ( 1/sum(market.gamma[sort_order[1:j-1]]) -
                      1/sum(market.gamma[sort_order[1:j]]) )
    end
    
    return A, sort_order
end


"""
    tatonnement(market)

Simultaneous tatonnement process for the given market.
"""
function tatonnement(market::Market; maxit::Int=50, p0=nothing,
                     rate=.5, damping=.001)
    if p0 === nothing 
        # Random cutoffs
        p = rand(length(market))
    else
        p = p0
    end
    
    res = Vector[]
    
    for it in 1:maxit 
        push!(res, p)
        p = max.(0, p + rate * (demand(market, p) - market.capacities) / it^damping)
    end
    
    return res
end


"""
    equilibrium(market)
    
Uses the gamma-over-capacity property to find the equilibrium.
"""
function equilibrium(market::Market)
    A, sort_order = demandmatrix(market, market.gamma ./ market.capacities)

    p_heuristic = zeros(length(market))
    p_heuristic[sort_order] = max.(0, A\(market.capacities - market.gamma/sum(market.gamma))[sort_order])
    
    return p_heuristic
end


"""
    makepreflists(market; n_students)

Generates preference lists for schools in the corresponding market.
"""
function makepreflists(market::Market; n_students::Int=100)
    caps = round.(Int, n_students * market.capacities)
    
    dist = Gumbel()
    
    n_schools = length(market)
    
    studentprefs = zeros(Int, n_schools, n_students)
    for i in 1:n_students
        ratings = market.qualities + rand(dist, n_schools)
        studentprefs[:, i] = invperm(sortperm(ratings, rev=true))
    end
    
    scores = rand(n_students)
    schoolprefs = repeat(invperm(sortperm(scores, rev=true)), 1, n_schools) 
    
    return studentprefs, schoolprefs, caps, scores
end

include("SchoolIncentives.jl")

end
