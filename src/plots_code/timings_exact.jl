using DataFrames, Distributions, BenchmarkTools, LaTeXStrings
setprecision(BigFloat, 128)
Ls = [L for L in 50:10:100]
Ns = [n for n in 10:20:100]

#warmup 
exact_randbow_clean(100, 1; save = false)

data = DataFrame(L = Int64[], N = Int64[], time_secs = Float64[])
for trio in Iterators.product(Ls, Ns)
    L, n_trials = trio
    tick = time()
    exact_randbow_clean(L, n_trials)
    tock = time()
    duration = tock - tick
    push!(data, (L, n_trials, duration))
end

