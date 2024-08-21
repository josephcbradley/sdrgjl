println("Randbow central negativity with OBC, varied r, XXY chain")
Δ = BigFloat(1.0)

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

plt = plot(dpi = 600, xlabel = L"r", ylabel = L"\mathcal{E}_{A_1:A_2}", 
    legend = :topright)

#params ordered δ, h
params = [(0.5, 1.5), (1.0, 3.0), (0.5, 2.), (1., 4.), (0.5, 3.5), (1., 7.)]

for ps in reverse(params)

    delta, h = ps
    randbow_distribution = RandbowDistribution{BigFloat}(delta, h)

    results = disorder_negativity_r_average(sdrg_params, 
        n_trials, randbow_distribution, l, min_window_size = min_window_size, 
        window_increment = window_increment, window_ceiling = window_ceiling,
        )

    #subset!(results, :negativity => x -> x .> 0)

    @df results scatter!(:l, :negativity, label = L"\delta = %$delta, h = %$h")
end

plt

saveplot(plt, img_path, "randbow_central_neg_obc_XXY_vary_r")