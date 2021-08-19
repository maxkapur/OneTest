cd ..
for i in 1 2 4 8
do
    julia --project --threads=$i ./BestResponseParallelization/BRP.jl
done

julia --project ./BestResponseParallelization/plot.jl
