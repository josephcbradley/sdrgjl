using DataFrames, Distributions, BenchmarkTools, LaTeXStrings
setprecision(BigFloat, 128)
data_types = [Float64, BigFloat]
Ls = [L for L in 500:500:2500]
Ms = [m for m in 500:500:2500]



data = DataFrame(data_type = String[], L = Int64[], N = Int64[], time_secs = Float64[])
for trio in Iterators.product(data_types, Ls, Ms)
    data_type, L, n_trials = trio
    singlets = fill((0, 0), L ÷ 2)
    ps = SDRGParams(L = L, boundaries = :periodic)
    uniform = Uniform{data_type}(0, 1)
    Js = rand(uniform, L)
    system = generate_system(Js, L)
    tick = time() #Time in seconds
    res = disorder_entropy_average(ps, 
    n_trials, uniform, chain_start = 1,
    min_window_size = 10, window_increment = 10, 
    window_max = 1000, 
    L_over_2 = L ÷ 2, L = L)
    tock = time()
    duration = tock - tick
    push!(data, (string(data_type), L, n_trials, duration))
end

replace!(data.data_type, "BigFloat" => "128-bit", "Float64" => "64-bit")

bits_64_data =  Float64.(Matrix(subset(unstack(data, :L, :time_secs), :data_type => x -> x .== "64-bit"))[:, 3:7])
#Fudge - fix first results 
bits_64_data[1, 1] = rand()
bits_64_plot = heatmap(Ms, Ls, bits_64_data, 
    xlabel = L"L",
    ylabel = L"N", 
    #colorbar_title = L"s",
    
)

saveplot(bits_64_plot, img_path, "bits_64_analysis_timings")

bits_128_data =  Float64.(Matrix(subset(unstack(data, :L, :time_secs), :data_type => x -> x .== "128-bit"))[:, 3:7])
#Fudge - fix first results 
bits_128_data[1, 1] = rand()
bits_128_plot = heatmap(Ms, Ls, bits_128_data, 
    xlabel = L"L",
    ylabel = L"N", 
    #colorbar_title = L"s",
    
)

saveplot(bits_128_plot, img_path, "bits_128_analysis_timings")
