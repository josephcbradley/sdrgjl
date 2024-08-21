println("Powerbow central negativity with OBC, varied r, XX chain")
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
    Δ = Δ)

plt = plot(dpi = 600, xlabel = L"r", ylabel = L"\mathcal{E}",
    legend = :topright)

#params ordered δ, alpha   
params = [(0.5, 1.5), (1.0, 3.0), (0.5, 2.), (1., 4.), (0.5, 3.5), (1., 7.)]

for ps in reverse(params)

    delta, h = ps
    powerbow_distribution = PowerbowDistribution{BigFloat}(delta, h)

    results = disorder_negativity_r_average(sdrg_params, 
        n_trials, powerbow_distribution, l, min_window_size = min_window_size, 
        window_increment = window_increment, window_ceiling = window_ceiling,
        )

    #subset!(results, :negativity => x -> x .> 0)

    @df results scatter!(:l, :negativity, label = L"\delta = %$delta, \alpha = %$h")
end

plt

saveplot(plt, img_path, "powerbow_central_neg_obc_XX_vary_r")