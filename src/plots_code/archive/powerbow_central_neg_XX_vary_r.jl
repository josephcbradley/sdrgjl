println("Randbow central negativity with OBC, varied r, XX chain")
Δ = BigFloat(0.0)

min_window_size = 2
window_increment = 2
window_ceiling = 100
l = 100
L = 1000
sdrg_params = SDRGParams(L = L,
    plots = false, 
    boundaries = :open, 
    verbose = false, 
    Δ = 0.0)

plt = plot(dpi = 600, xlabel = L"r", ylabel = "\mathcal{E}", yscale = :log10)

#params ordered δ, h
params = [(0.5, 1.5), (1.0, 3.0), (1.5, 4.5)]

for ps in params

    delta, alpha = ps
    randbow_distribution = powerbowDistribution{BigFloat}(delta, alpha)

    results = disorder_negativity_r_average(sdrg_params, 
        n_trials, randbow_distribution, l, min_window_size = min_window_size, 
        window_increment = window_increment, window_ceiling = window_ceiling,
        start_pos = (L ÷ 2) - min_window_size)

    subset!(results, :negativity => x -> x .> 0)

    @df  results scatter(:l, :negativity, label = L"\delta = %$delta, \alpha = %$alpha")
end

saveplot(plt, img_path, "randbow_central_neg_obc_XX_vary_r")