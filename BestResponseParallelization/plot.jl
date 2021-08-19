using Plots, CSV, DataFrames, Optim

function fit_powerlaw_xyspace(x, y; b0=4)
    n = length(x)
    
    # z = [a, b]
    f(z) = sum( (z[1] * x[i] ^ z[2] - y[i])^2 for i in 1:n)
    
    a, b = optimize(f, [1e-5, b0], autodiff=:forward).minimizer
   
    return a, b
end

pl = plot(xlabel="number of schools",
          ylabel="best-response iteration time",
          legend=:topleft,
          xlim=(15, 155))

fnames = readdir("./BestResponseParallelization/times/")
colors = theme_palette(:auto)



for (i, fname) in enumerate(fnames)
    df = DataFrame(CSV.File("./BestResponseParallelization/times/"*fname))

    lab = fname[begin:end-11]

    if lab=="1"
        lab *= " thread"
    else
        lab *= " threads"
    end

    # Fit a power law: time = exp(a) * C^p
    # a, p = hcat(ones(length(df.C)), log.(df.C)) \ log.(df.time)
    # a = exp(a)
    a, p = fit_powerlaw_xyspace(df.C, df.time)

    scatter!(pl, df.C, df.time, label=lab, c=colors[i], msw=0, alpha=0.4, ms=3)
    plot!(pl, x->a*x^p, c=colors[i], label="")
end

savefig(pl, "./BestResponseParallelization/results.pdf")
