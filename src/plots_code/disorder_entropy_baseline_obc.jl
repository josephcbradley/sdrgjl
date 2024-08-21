# Entanglement Entropy spectrum w/ OBC? 
params_1000_open = SDRGParams(L = 1000, plots = false, 
    boundaries = :open, img_path = img_path, 
    data_path = data_path, Δ = 1.0)
results_1000_open = disorder_entropy_average(params_1000_open, 
    n_trials, uniform_distribution, chain_start = 1,
    min_window_size = 10, window_increment = 10, 
    window_max = 1000, 
    L_over_2 = 500, L = 1000)
adjust_l_entropy_data!(results_1000_open, params_1000_open.L)

params_2000_open = SDRGParams(L = 2000, plots = false, 
    boundaries = :open, img_path = img_path, data_path = data_path,
    Δ = 1.0)

results_2000_open = disorder_entropy_average(params_2000_open, n_trials, 
    uniform_distribution, min_window_size = 10, window_increment = 10, 
    L_over_2 = 1000, L = 2000, chain_start = 1, window_max = 2000)
adjust_l_entropy_data!(results_2000_open, params_2000_open.L)
#Create logged l as well
results_2000_open.l_adj_log10 .= log.(Ref(10), results_2000_open.l_adj)

# dodgy numbers coming in - filter appropriatley 
filter!(:l_adj_log10 => x -> x > zero(x), results_2000_open)


obc_results_plt = @df results_2000_open scatter(:l_adj, :entropy, markercolor = :green,
    markershape = :square, alpha = 0.8, xscale = :log10, markersize = 5.0,
    label = L"L = 2000", legend = :bottomright, markerstrokewidth = 2, 
    markerstrokecolor = :green, dpi = 600,
    ylims = (0.7, 1.3),
    xlims = (10, 2000),
    xlabel = L"\frac{L}{\pi} Y\left(\frac{\pi l}{L}\right)", ylabel = L"S_A",
    clip_mode = "individual")

@df results_1000_open scatter!(obc_results_plt, :l_adj, :entropy, color = :blue,
    label = L"L = 1000", alpha = 0.8, markerstrokecolor = :blue, markersize = 5.0,
    clip_mode = "individual")


ols = lm(@formula(entropy ~ l_adj_log10), results_2000_open)

plot!(obc_results_plt, results_2000_open.l_adj, predict(ols), label = "SDRG",
    colour = :red, linestyle = :dashdot)

saveplot(obc_results_plt, img_path, "entropy_baseline_obc")

obc_results_plt = @df results_2000_open scatter(:l, :entropy, markercolor = :green,
    markershape = :square, alpha = 0.8, xscale = :log10, markersize = 5.0,
    label = L"L = 2000", legend = :bottomright, markerstrokewidth = 2, 
    markerstrokecolor = :green, dpi = 600,
    ylims = (0.7, 1.3),
    xlims = (10, 2000),
    xlabel = L"l", ylabel = L"S_A")

@df results_1000_open scatter!(obc_results_plt, :l, :entropy, color = :blue,
    label = L"L = 1000", alpha = 0.8, markerstrokecolor = :blue, markersize = 5.0)

saveplot(obc_results_plt, img_path, "entropy_baseline_obc_raw")