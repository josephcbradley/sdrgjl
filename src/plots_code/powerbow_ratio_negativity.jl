L = 240
Δ = BigFloat(0)
ps = SDRGParams(L = L, plots = false, 
    boundaries = :open, img_path = img_path, 
    data_path = data_path, Δ = Δ)
L_over_2 = L ÷ 2
hds = convert.(Tuple{BigFloat, BigFloat}, [(7, 1), (14, 2), (3.5, 0.5), (4, 1), (8, 2), (2, .5), (0.5, 1), (0.25, 0.5)])
min_window_size = 1
window_increment = 1
window_adjust = -1
r = 0
start_pos = (L - r) ÷ 2
fig_6 = plot(dpi = 600, xlabel = "l", ylabel = "S",
    #xlims = (2, 100 + 10), 
    #xscale = :log10, 
    #yscale = :log10, 
    legend = :topleft)

for pair in reverse(hds)
    α, δ = pair

    powerbow = PowerbowDistribution(δ, α)

    disd = disorder_negativity_average(ps, n_trials, powerbow, r,
        min_window_size = min_window_size,
        window_increment = window_increment,
        window_adjust = window_adjust,
        #L_over_2 = L_over_2, 
        #window_max = L_over_2,
        L = L,
        start_pos = start_pos)

    @df disd scatter!(fig_6, :l, :negativity, label = "α = $α, δ = $δ")

    # Fit square root 
    optim_result = optimize(x -> sqrt_loss(x, disd.l, disd.negativity), ones(2))
    a, b = Optim.minimizer(optim_result)
    plot!(fig_6, x -> a + (b*sqrt(x)), 10:0.1:100, label = nothing,
    colour = :red, linestyle = :dashdot)

end

fig_6

saveplot(fig_6, img_path, "randbow_ratio_entropy")