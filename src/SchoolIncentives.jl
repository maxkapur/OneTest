const phi = Float16((-1 + √5)/2)
const omp = Float16((3 - √5)/2) # One minus phi


"""
    utility(market, sigma, p)

Gives the value of each school's objective function ``u_c(p_c) = \\log D_c + σ \\log p_c``.
"""
function utility(market::Market{T},
                 sigma::Union{<:Integer, <:AbstractFloat, Vector{<:AbstractFloat}},
                 p::Vector{<:AbstractFloat})::Vector{T} where T<:AbstractFloat
    return log.(demand(market, p)) .+ sigma .* log.(p)
end


"""
   incentivegradient(market, sigma, cutoffs)

When each school's objective function is ``u_c(p_c) = \\log D_c + σ \\log p_c``,
gives the gradient of each school's objective with respect to its own cutoff.
"""
function incentivegradient(mkt::Market,
                           sigma::Union{<:Integer, <:AbstractFloat, Vector{<:AbstractFloat}},
                           p::Vector{<:AbstractFloat})
    A, sort_order = demandmatrix(mkt, p)
    
    grad = zeros(length(mkt))
    grad[sort_order] = diag(A) ./ demand(mkt, p)[sort_order] .+ sigma[sort_order] ./ p[sort_order]
    
    return grad
end


"""
    localequilibriumsearch(mkt, sigma; maxit, p0, rate, damping)

Search for a local equilibrium by having each school follow its incentive gradient
according to decreasing sequence of step sizes.
"""
function localequilibriumsearch(mkt::Market, sigma::Union{<:Integer, <:AbstractFloat, Vector{<:AbstractFloat}}; maxit::Int=50, p0=nothing, rate=.005, damping=.2)
    if p0===nothing
        p = rand(length(mkt))
    else
        p = p0
    end
    
    res = Vector[]
    clamp = 1e-4
    
    for it in 1:maxit 
        push!(res, p)
        grad = incentivegradient(mkt, sigma, p)
        p = min.(1 - clamp, max.(clamp, p + rate * grad / it^damping))
    end
    
    return res
end


"""
    goldensection(f, lb, ub, nit)

Computes the maximizer of `f` on the interval `(lb, ub)` by the method of golden
sections. Uses `nit` iterations. Currently unused.
"""
function goldensection(f::Function,
                       lb::T,
                       ub::T;
                       nit=25, verbose=false)::T where T<:AbstractFloat
    verbose && print("j = ")
    for j in 1:nit
        verbose && print("$j ")
        a = lb + omp*(ub - lb)
        b = lb + phi*(ub - lb)

        if f(a) < f(b)
            lb = a
        else
            ub = b
        end
    end
    
    return (lb + ub) / 2
end


"""
    brforc(c, sort_order, market, sigma, p; heuristic=false)

Computes the best response for school `c`. If `heuristic=true`, uses the
"cloud hopping" heuristic, which reduces the computation time if the local optima
form a discrete quasiconvex function. 
"""
@inline function brforc(c::Int,
                        sort_order::Vector{Int},
                        market::Market{T},
                        sigma::T,
                        p::Vector{T};
                        heuristic=false::Bool
                        )::T where T<:AbstractFloat
    
    p = vcat(p, 1.)
    C = length(market)
    
    global_sort_order = copy(sort_order)
    local_sort_order = vcat(global_sort_order, C+1)
    
    setdiff!(global_sort_order, c)
    
    gamma_cumsum = cumsum(market.gamma[global_sort_order])
    
    if heuristic
        @inline function fetch_̄p(i::Int)::Tuple{T, T}
            if i > 1
                local_sort_order[begin:i-1] = global_sort_order[begin:i-1]
            end
            local_sort_order[i] = c
            local_sort_order[i+1:end-1] = global_sort_order[i:end]

            if i==C
                density_cumsum = 0
                p̄ = sigma / (1 + sigma) # * (p[C+1] - density_cumsum) = 1
            else
                density_cumsum = sum((p[local_sort_order[d+1]] - p[local_sort_order[d]]) /
                                     (market.gamma[c] + get(gamma_cumsum, d-1, 0))
                                     for d in i+1:C)

                p̄ = (sigma / (1 + sigma)) * 
                    (p[local_sort_order[i+1]] +
                        (market.gamma[c] + get(gamma_cumsum, i-1, 0)) * density_cumsum)
            end

            # clamp
            lb = i==1 ? 0 : p[local_sort_order[i-1]]
            ub = p[local_sort_order[i+1]]
            p̄ = max(lb, min(ub, p̄))
            D̄ = density_cumsum + (p[local_sort_order[i+1]] - p̄) / 
                                 sum(market.gamma[local_sort_order[1:i]])
            # gamma[c] * D̄ = demand for c

            Ū = D̄ * p̄^sigma
            
            return p̄, Ū
        end
        
        
        p̄_vec = zeros(T, C+1)
        Ū_vec = fill(T(Inf), C+1)
        Ū_vec[C+1] = 0

        h_lb, h_ub = 1, C
        h = 0

        it = 0
        while h_ub - h_lb ≥ 1 && it < 1000
            it += 1
            h = h_ub + (h_lb - h_ub) ÷ 2
                
            if Ū_vec[h]==Inf
                p̄_vec[h], Ū_vec[h] = fetch_̄p(h)
            end

            if Ū_vec[h+1]==Inf
                p̄_vec[h+1], Ū_vec[h+1] = fetch_̄p(h+1)
            end

            if Ū_vec[h+1] - Ū_vec[h] ≤ 0
                if h == 1
                    # Descending at first point: Optimal
                    h_lb = h
                else
                    # Check the previous index
                    if Ū_vec[h-1]==Inf
                        p̄_vec[h-1], Ū_vec[h-1] = fetch_̄p(h-1)
                    end 
                    if Ū_vec[h] - Ū_vec[h-1] > 0
                        # Found low-high-low pattern: Optimal
                        h_lb = h
                    end
                end
                h_ub = h
            else
                h_lb = h
            end    
        end
        
        return p̄_vec[h]

    else   # Search all indices
        p_star = T(0)
        U_star = T(0)
        
        for i in 1:C
            if i > 1
                local_sort_order[i-1] = global_sort_order[i-1]
            end
            local_sort_order[i] = c
            local_sort_order[i+1:end-1] = global_sort_order[i:end]

            if i==C
                density_cumsum = 0
                p̄ = sigma / (1 + sigma) # * (p[C+1] - density_cumsum) = 1
            else
                density_cumsum = sum((p[local_sort_order[d+1]] - p[local_sort_order[d]]) /
                                     (market.gamma[c] + get(gamma_cumsum, d-1, 0))
                                     for d in i+1:C)

                p̄ = (sigma / (1 + sigma)) * 
                    (p[local_sort_order[i+1]] +
                        (market.gamma[c] + get(gamma_cumsum, i-1, 0)) * density_cumsum)
            end

            # clamp
            lb = i==1 ? 0 : p[local_sort_order[i-1]]
            ub = p[local_sort_order[i+1]]
            p̄ = max(lb, min(ub, p̄))
            D̄ = density_cumsum + (p[local_sort_order[i+1]] - p̄) / 
                                 sum(market.gamma[local_sort_order[1:i]])
            # gamma[c] * D̄ = demand for c

            Ū = D̄ * p̄^sigma

            if Ū > U_star
                p_star, U_star = p̄, Ū
            end
        end

        return p_star
    end
end


"""
    bestresponse_it(market, sigma, p; verbose=false, heuristic=false)

Computes one iteration of best_response dynamics, wherein each school
updates its cutoff to the value that maximizes its utility on the assumption
that other schools' cutoffs remain fixed. If `heuristic=true`, use the
"cloud hopping" heuristic, which reduces the computation time if the local optima form a discrete quasiconvex function (as is usually true when the number)
of schools is larger than about 10). 
"""
function bestresponse_it(market::Market{T},
                         sigma::Array{T, 1},
                         p=nothing::Union{Nothing, Array{T, 1}};
                         verbose=false::Bool,
                         heuristic=false::Bool
                         )::Vector{T} where T<:AbstractFloat
    
    C = length(market)
    global_sort_order = sortperm(p)
    
    p_star = zeros(T, C)
    
    verbose && print("\n  c = ")
    
    for c in 1:C
        verbose && c%10 == 1 && print("$c ")
        p_star[c] = brforc(c, global_sort_order, market, sigma[c], p; heuristic=heuristic)
    end
    
    return p_star
end


"""
    bestresponse(market, sigma, p_curr; nit=10)

Compute iterates of the market when each school adjusts its cutoff to maximize its
utility, where each school's utility function is ``u_c(p_c) = \\log D_c + σ \\log p_c``.
"""
function bestresponse(market::Market{T},
                      sigma::Array{T, 1},
                      p_curr=nothing::Union{Nothing, Array{T, 1}};
                      nit=10::Int
                      )::Vector{Vector{T}} where T<:AbstractFloat
    if p_curr===nothing
        p_curr = rand(length(mkt))
    end
    
    res = Vector[p_curr]
    
    for j in 1:nit
        p_curr = bestresponse_it(market, sigma, p_curr)
        push!(res, copy(p_curr))
    end
    
    return res
end


"""
    secantroot(f, x1, x2, nit)

Use the secant method to find roots of `f` given starting points x1 and x2.
"""
function secantroot(f::Function,
                    x1::T,
                    x2::T;
                    maxit=25, epsilon=Float16(1e-6), verbose=false
                    )::Tuple{T, T} where T<:AbstractFloat
    i = 0
    y1, y2 = f(x1), f(x2)
    
    while i < maxit && abs(y2) > epsilon && !(x1 == x2)
        i += 1
        x3 = x2
        
        # Prevent NaN if we hit the zero exactly
        if (dy = y2 - y1) != 0
            x3 = max(0 + epsilon, x1 - y1 * (x2 - x1) / dy)
        end
        
        x1, x2 = x2, x3
        y1, y2 = y2, f(x3)
    end
    
    verbose && println("  Iterations in secant: $i")
    
    return x2, y2
end


"""
    bisection(f, x1, x2, nit)

Use the bisection method to find roots of `f` given starting points `x1` and `x2`. `x1` must be
a lower bound; `x2` will be doubled until an upper bound is discovered.
"""
function bisection(f::Function,
                    x1::T,
                    x2::T;
                    maxit=25,
                    epsilon=Float16(1e-6),
                    verbose=false
                    )::Tuple{T, T} where T<:AbstractFloat

    max_ub_search = 25
    y1, y2 = f(x1), f(x2)
    
    h = 0
    while y2 < 0 && h < max_ub_search
        h += 1
        x2 *= 2
        y2 = f(x2)
    end
    
    h == max_ub_search ? @warn("Failed to find an upper bound after $h doublings") : nothing
    
    x3, y3 = x2, y2
    
    i = 0
    while i < maxit && abs(y2) > epsilon && !(x1 == x2)
        i += 1
        
        x3 = (x1 + x2)/2
        y3 = f(x3)
        
        if y3 < 0
            # New LB; keep UB. 
            x1, x2 = x3, x2
            y1, y2 = y3, y2
        else
            # New UB; keep LB. 
            x1, x2 = x1, x3
            y1, y2 = y1, y3
        end
    end
    
    verbose && println("  Iterations in bisection: $i")
    
    return x3, y3
end


"""
    sigmainvopt(market, p; x1, x2, maxit=25, epsilon=1e-5, verbose=false, heuristic=false)

Estimate the selectivity parameter ``σ`` for each school such that the given cutoffs are 
an equilibrium, when school's objective functions are ``u_c(p_c) = \\log D_c + σ \\log p_c``.
Returns the vector of ``σ``-values and an error vector.

If `heuristic=true`, uses the
"cloud hopping" heuristic (see `?bestresponse_it`). 
"""
function sigmainvopt(market::Market{T},
                     p::Array{T, 1};
                     x1=nothing::Union{Nothing, Array{T, 1}},
                     x2=nothing::Union{Nothing, Array{T, 1}}, 
                     maxit=10::Int,
                     epsilon=Float16(1e-5),
                     verbose=false::Bool,
                     heuristic=false::Bool
                     )::Tuple{Vector{T}, Vector{T}} where T<:AbstractFloat
    
    if x1 === nothing
        x1 = fill(T(0), length(market))
    end
    
    if x2 === nothing
        x2 = fill(T(2e-1), length(market))
    end
    
    global_sort_order = sortperm(p)
    
    res = zeros(T, length(market))
    err = zeros(T, length(market))
    
    for (c, pc) in enumerate(p)
        verbose && @show c
        if pc == 0
            res[c] = 0
        else
            function f(x::T)::T
                return brforc(c, global_sort_order, market, x, p,
                              heuristic=heuristic) - p[c]
            end
            
            res[c], err[c] = bisection(f, x1[c], x2[c];
                                       maxit=maxit,
                                       epsilon=epsilon,
                                       verbose=verbose)
        end
    end
    
    return res, err
end
