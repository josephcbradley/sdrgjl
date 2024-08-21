println("Running baseline entanglement negativity plot for r = 0, OBC, with uniform distribution.")
params_1000_open = SDRGParams(L = 1000, plots = false, 
    boundaries = :open, img_path = img_path, 
    data_path = data_path, Δ = 1.0)

results_1000 = disorder_negativity_average(params_1000_open, n_trials, uniform_distribution, 0,
    start_pos = 1, min_window_size = 10, window_increment = 10,
    L = 1000)
add_shifted_negativity!(results_1000, params_1000_open.L)
#also do l/L
results_1000.l_adj = results_1000.l ./ Ref(1000)

params_2000_open = SDRGParams(L = 2000, plots = false, 
    boundaries = :open, img_path = img_path, 
    data_path = data_path, Δ = 1.0)

results_2000 = disorder_negativity_average(params_2000_open, n_trials, uniform_distribution, 0,
    start_pos = 1, min_window_size = 10, window_increment = 10,
    L = 2000)
add_shifted_negativity!(results_2000, params_2000_open.L)
results_2000.l_adj = results_2000.l ./ Ref(2000)


results_plt = @df results_2000 scatter(:l_adj, :shifted_negativity, markercolor = :green,
    markershape = :square, alpha = 0.8, markersize = 5.0,
    label = L"L = 2000", legend = :bottomright, markerstrokewidth = 2, 
    markerstrokecolor = :green, dpi = 600, 
    #ylims = (1.0, 2.25),
    xlabel = L"l/L", ylabel = L"\mathcal{E}^s_{A_1:A_2}",
    clip_mode = "individual")

@df results_1000 scatter!(results_plt, :l_adj, :shifted_negativity, color = :blue,
label = L"L = 1000", alpha = 0.8, markerstrokecolor = :blue, markersize = 5.0,
clip_mode = "individual")

optim_result = optimize(x -> infinite_chain_negativity_loss(x, results_2000.l, results_2000.negativity), ones(1))
k = Optim.minimizer(optim_result)[1] - 0.03 #fudge factor
plot!(l -> infinite_chain_negativity(l, k), 0.0:0.01:0.5, label = "SDRG",
colour = :red, linestyle = :dashdot)

saveplot(results_plt, img_path, "baseline_negativity_plot_obc")

CSV.write(data_path * "baseline_negativity_data_obc" * "_1000.csv", results_1000)   
CSV.write(data_path * "baseline_negativity_data_pbc" * "_2000.csv", results_2000)