# Build baseline entanglement entropy spectrum plots
println("Running baseline entropy spectrum plots with uniform distribution.")

params_1000_periodic = SDRGParams(L = 1000, plots = false, 
    boundaries = :periodic, img_path = img_path, 
    data_path = data_path, Δ = 1.0)

results_1000 = disorder_entropy_average(params_1000_periodic, n_trials, 
    uniform_distribution,
    min_window_size = 10, window_increment = 10, chain_start = 1,
    L_over_2 = 500, L = 1000)
adjust_l_entropy_data!(results_1000, params_1000_periodic.L)

params_2000_periodic = SDRGParams(L = 2000, plots = false, 
    boundaries = :periodic, img_path = img_path, 
    data_path = data_path, Δ = 1.0)

    results_2000 = disorder_entropy_average(params_2000_periodic, n_trials, uniform_distribution,
    min_window_size = 10, window_increment = 10, chain_start = 1,
    L_over_2 = 1000, L = 2000)
adjust_l_entropy_data!(results_2000, params_2000_periodic.L)
#Create logged l as well
results_2000.l_adj_log10 .= log.(Ref(10), results_2000.l_adj)

results_plt = @df results_2000 scatter(:l_adj, :entropy, markercolor = :green,
    markershape = :square, alpha = 0.8, xscale = :log10, markersize = 5.0,
    label = L"L = 2000", legend = :bottomright, markerstrokewidth = 2, 
    markerstrokecolor = :green, dpi = 600,
    # ylims = (1.0, 2.25), 
    xlabel = L"\frac{L}{\pi} Y\left(\frac{\pi l}{L}\right)", ylabel = L"S_A",
    clip_mode = "individual")

@df results_1000 scatter!(results_plt, :l_adj, :entropy, color = :blue,
    label = L"L = 1000",
     alpha = 0.8, markerstrokecolor = :blue, markersize = 5.0, clip_mode = "individual")

#Add regression line
ols = lm(@formula(entropy ~ l_adj_log10), results_2000)

plot!(results_plt, results_2000.l_adj, predict(ols), label = "SDRG",
    colour = :red, linestyle = :dashdot)

saveplot(results_plt, img_path, "entropy_baseline_pbc")

results_plt = @df results_2000 scatter(:l, :entropy, markercolor = :green,
    markershape = :square, alpha = 0.8, xscale = :log10, markersize = 5.0,
    label = L"L = 2000", legend = :bottomright, markerstrokewidth = 2, 
    markerstrokecolor = :green, dpi = 600,
    # ylims = (1.0, 2.25),
    xlabel = L"l", ylabel = L"S_A")

@df results_1000 scatter!(results_plt, :l, :entropy, color = :blue,
    label = L"L = 1000", alpha = 0.8, markerstrokecolor = :blue, markersize = 5.0)



saveplot(results_plt, img_path, "entropy_baseline_pbc_raw")



#= CSV.write(data_path * "entropy_baseline_pbc" * "_1000.csv", results_1000)   
CSV.write(data_path * "entropy_baseline_pbc" * "_2000.csv", results_2000) =#