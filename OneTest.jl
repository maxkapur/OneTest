using LinearAlgebra
import Base:length

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
    Market(qualities, capacities)

Instantiate a school-choice market. `gamma` is normalized to sum to one.
"""
function Market(qualities::Array{T, 1}, capacities::Array{T, 1}) where {T<:AbstractFloat}
    length(qualities) == length(capacities) || return throw(DimensionMismatch)
    gamma = exp.(qualities)
    gamma /= sum(gamma)
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
    appeal(market, cutoffs)

Return appeal of entring class at each school given a set of cutoffs, using
multinomial logit choice model with one student profile and a single test score.
"""
function appeal(market      ::Market,
                cutoffs     ::AbstractArray{<:AbstractFloat, 1};
                )::AbstractArray{<:AbstractFloat, 1}
    
    m = length(market)
    demands = zeros(m)

    sort_order = sortperm(cutoffs)
    cutoffs[sort_order]

    γ = exp.(market.qualities)
    appeals = zeros(m)

    diff_of_squares = diff([cutoffs[sort_order]; 1] .^ 2)

    for c in 1:m, d in c:m     # For each score threshold
        appeals[sort_order[c]] += diff_of_squares[d] *
                                  γ[sort_order[c]] / sum(γ[sort_order[1:d]])
    end

    return appeals / 2
end


"""
    demandmatrix(market, cutoffs)

Returns the matrices `A` and the permutation `sort_order` such that

````
demand(market, p)[sort_order] == A * cutoffs[sort_order] + market.gamma[sort_order]
````

and

````
appeal(market, p)[sort_order] == (A * cutoffs[sort_order] .^2 + market.gamma[sort_order])/2
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
    quickeq(market)

Tatonnement procedure for finding equibilibrium.
"""
function equilibrium(market::Market; maxit::Int=500, tol=1e-8, p0=nothing)
    if p0 === nothing || p0 == :heuristic
        # Heuristic cutoffs
        p = heuristiceq(market)
    elseif p0 == :random
        # Random cutoffs
        p = rand(length(market))
    else
        p = p0
    end
    
    nit = 0
    
    while true
        nit += 1
        D = demand(market, p)
        
        if all(D .< market.capacities .+ tol) || nit==maxit
            break
        end
#         Generic decreasing step sequence
        p = max.(0, p + .5 * (D - market.capacities) / nit^.001)
        
#         Step sequence based on monopoly assumption
#         Seems clever but can cause cycling in practice
#         p = max.(0, p + (1 .- market.capacities ./ D) .* (1 .- p))
        
#         Normalize to satisfy D = Ap + γ
        p = heuristiceq(market; orderby=p)
    end
    
    return p, nit
end


"""
    heuristiceq(market; orderby)
    
Uses the gamma-minus-capacity heuristic to find an approximate equilibrium, or tries to order the optimal cutoffs by another argument supplied. 
"""
function heuristiceq(market::Market; orderby=nothing)
    if orderby === nothing
        orderby = market.gamma - market.capacities
    end
    
    # Use the heuristic as the cutoffs arg, since demand() looks only at its order
    A, sort_order = demandmatrix(market, orderby)

    p_heuristic = zeros(length(market))
    p_heuristic[sort_order] = max.(0, A\(market.capacities - market.gamma)[sort_order])
    
    return p_heuristic
end