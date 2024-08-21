#Now for XXY 
hs = (BigFloat.([0.5, 3., 7., 10.]))
negativity_randbow_plot_vec = Plots.Plot{Plots.GRBackend}[]
for h in reverse(hs)
    println("Running entanglement negativity plot for r = 0 with h = $h, δ = 1.0")
    randbow = RandbowDistribution(BigFloat(1.0), BigFloat(h))
    randbow_params = SDRGParams(L = 1000, boundaries = :open, Δ = 1.0)
    randbow_results_1000 = disorder_negativity_average(randbow_params, ntrials, randbow, 0,
        start_pos = 1, min_window_size = 10, window_increment = 10,
        L = randbow_params.L)
    add_shifted_negativity!(randbow_results_1000, randbow_params.L)
    #also do l/L
    randbow_results_1000.l_adj = randbow_results_1000.l ./ Ref(randbow_params.L)

    randbow_params = SDRGParams(L = 2000, boundaries = :open, Δ = 1.0)
    randbow_results_2000 = disorder_negativity_average(randbow_params, ntrials, randbow, 0,
        start_pos = 1, min_window_size = 10, window_increment = 10,
        L = randbow_params.L)
    #add_shifted_negativity!(randbow_results_2000, randbow_params.L)
    randbow_results_2000.l_adj = randbow_results_2000.l ./ Ref(2000)


    results_plt = @df randbow_results_2000 scatter(:l_adj, :negativity, markercolor = :green,
        markershape = :square, alpha = 0.8, markersize = 4.0,
        label = L"L = 2000", legend = :topleft, markerstrokewidth = 1, 
        markerstrokecolor = :green, dpi = 600, title = "h = $h, δ = 1.0",
        #ylims = (1.0, 2.25),
        xlabel = L"l/L", ylabel = L"\mathcal{E}_{A_1:A_2}")

    @df randbow_results_1000 scatter!(results_plt, :l_adj, :negativity, color = :blue,
    label = L"L = 1000", alpha = 0.8, markerstrokecolor = :blue, markersize = 4.0)

    push!(negativity_randbow_plot_vec, results_plt)
end

fixed_d_vary_h_negativity_delta_1 = plot(negativity_randbow_plot_vec..., layout = (2, 2), dpi = 600)

savefig(fixed_d_vary_h_negativity_delta_1, img_path * "fixed_d_vary_h_negativity_delta_1" * ".png")