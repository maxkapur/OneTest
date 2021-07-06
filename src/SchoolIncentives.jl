const phi = (-1 + √5)/2
const omp = (3 - √5)/2 # One minus phi


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
function incentivegradient(mkt::Market, sigma::Union{<:Integer, <:AbstractFloat, Vector{<:AbstractFloat}},
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
sections. Uses `nit` iterations.
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
    bestresponse(market, sigma, p_curr; nit, golden_nit)

Compute iterates of the market when each school adjusts its cutoff to maximize its
utility, where each school's utility function is ``u_c(p_c) = \\log D_c + σ \\log p_c``.
Set `golden_nit=0` to maximize over the set of cutoffs instead of using line search.
"""
function bestresponse(market::Market{T},
                      sigma::Array{T, 1},
                      p_curr=nothing::Union{Nothing, Array{T, 1}};
                      nit=10::Int, golden_nit=5::Int, verbose=false
                      )::Vector{Vector{T}} where T<:AbstractFloat
    if p_curr===nothing
        p_curr = rand(length(mkt))
    end
    
    res = Vector[p_curr]
    
    for j in 1:nit
        p_next = zeros(length(p_curr))
        verbose && print("c = ")
        
        p_intervals = sort(p_curr)
        
        for c in 1:length(market)
            verbose && print("$c ")
            
            @inline function u(x::T)::T
                pc = p_curr[c]
                p_curr[c] = x
                out = utility(market, sigma, p_curr)[c]
                p_curr[c] = pc
                return out
            end
            
            if golden_nit == 0
                p_next[c] = argmax(u, p_intervals)
            else
                p_next[c] = argmax(u, goldensection(u, lb, ub, nit=golden_nit)
                                      for (lb, ub) in zip(vcat(0., p_intervals), vcat(p_intervals, 1.)))
            end
        end
        verbose && println()
        
        p_curr = copy(p_next)
        push!(res, p_curr)
    end
    
    return res
end


"""
    secantroot(f, x1, x2, nit)

Use the secant method to find roots of `f` given starting points x1 and x2.
Assume `f` is vector-valued where each entry of `f` is univariate in the corresponding `x[i]`. 
"""
function secantroot(f::Function,
                    x1::Array{T, 1},
                    x2::Array{T, 1};
                    nit=25, epsilon=1e-6, verbose=false)::Vector{T} where T<:AbstractFloat
    i = 0
    while i < nit && !(x1 == x2)
        i += 1
        verbose && @show i
        y1, y2 = f(x1), f(x2)
        x3 = copy(x2)
        for (i, dy) in enumerate(y2 .- y1)
            # Prevent NaN if we hit the zero exactly
            if dy != 0
                x3[i] = max.(0+epsilon, x1[i] - y1[i] * (x2[i] - x1[i]) / dy)
            end
        end
        x1, x2 = x2, x3
    end
    
    return (x1 + x2) ./ 2
end


"""
    sigmainvopt(market, p, secant_nit, golden_nit)

Estimate the selectivity parameter ``σ`` for each school such that the given cutoffs are 
an equilibrium, when school's objective functions are ``u_c(p_c) = \\log D_c + σ \\log p_c``.
"""
function sigmainvopt(market::Market{T}, p::Array{T, 1};
                     secant_nit=0::Int, golden_nit=0, verbose=false) where T<:AbstractFloat
    if secant_nit == 0
        secant_nit = length(market)
    end
    
    x1 = fill(T(1e-1), length(market))
    x2 = fill(T(1e0), length(market))
    
    @inline function f(x::Vector{T})::Vector{T}
        return p - bestresponse(market, x, p; nit=1, golden_nit=golden_nit, verbose=verbose)[end]
    end
    
    res = secantroot(f, x1, x2, nit=secant_nit, verbose=verbose)
    
    for (i, pc) in enumerate(p)
        if pc == 0
            res[i] = 0.
        end
    end
    
    return res
end
