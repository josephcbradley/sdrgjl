println("Randbow central negativity with OBC, XXY chain")
Δ = BigFloat(1.0)

min_window_size = 1
window_increment = 1
window_adjust = -1

randbow_distribution = RandbowDistribution{BigFloat}(1.0, 10.0)
params_1000_open_XX = SDRGParams(L = 1000, plots = false, 
    boundaries = :open, img_path = img_path, 
    data_path = data_path, Δ = Δ)

results_1000 = disorder_negativity_average(params_1000_open_XX, n_trials, 
    randbow_distribution, 0, start_pos = (params_1000_open_XX.L ÷ 2), 
    min_window_size = min_window_size, window_increment = window_increment,
    window_adjust = window_adjust,
    L = params_1000_open_XX.L)

add_shifted_negativity!(results_1000, params_1000_open_XX.L)
#also do l/L
results_1000.l_adj = results_1000.l ./ Ref(1000)

params_2000_open_XX = SDRGParams(L = 2000, plots = false, 
    boundaries = :open, img_path = img_path, 
    data_path = data_path, Δ = Δ)

results_2000 = disorder_negativity_average(params_2000_open_XX, n_trials, 
    randbow_distribution, 0, start_pos = (params_2000_open_XX.L ÷ 2), 
    min_window_size = min_window_size, window_increment = window_increment,
    window_adjust = window_adjust,
    L = params_2000_open_XX.L)
add_shifted_negativity!(results_2000, params_2000_open_XX.L)
results_2000.l_adj = results_2000.l ./ Ref(2000)

#Fudge - trim final results
filter!(:l => x -> x  < maximum(results_1000.l) - 2, results_1000) 
filter!(:l => x -> x  < maximum(results_2000.l) - 2, results_2000) 

results_plt = @df results_2000 scatter(:l_adj, :negativity, markercolor = :green,
    markershape = :square, alpha = 0.8, markersize = 5.0,
    label = L"L = 2000", legend = :bottomright, markerstrokewidth = 2, 
    markerstrokecolor = :green, dpi = 600, 
    #ylims = (1.0, 2.25),
    clip_mode = "individual",
    xlabel = L"l/L", ylabel = L"\mathcal{E}_{A_1:A_2}")

@df results_1000 scatter!(results_plt, :l_adj, :negativity, color = :blue,
    label = L"L = 1000", alpha = 0.8, markerstrokecolor = :blue,
    markersize = 5.0,
    clip_mode = "individual")


saveplot(results_plt, img_path, "randbow_central_neg_obc_XXY")
#=
CSV.write(data_path * "baseline_negativity_data_pbc" * "_1000.csv", results_1000)   
CSV.write(data_path * "baseline_negativity_data_pbc" * "_2000.csv", results_2000) =#