println("Randbow ratio entropy plot")
L = 240
ps = SDRGParams(L, false, :open, false, "", "", 0.0)
hds = convert.(Tuple{BigFloat, BigFloat}, [(7, 1), (14, 2), (3.5, 0.5), (4, 1), (8, 2), (2, .5), (0.5, 1), (0.25, 0.5)])
chain_start = L_over_2 + 1
min_window_size = 2
window_increment = 1
l_range = min_window_size:window_increment:L_over_2
l_length = length(l_range)

fig_6 = plot(dpi = 600, xlabel = "l", ylabel = "S",
    xlims = (2, 100 + 10), 
    xscale = :log10, 
    yscale = :log10, 
    #legend = :topleft # Legend seems to be obscuring all the data with the 
    # new font size, so let PGFPlots hang out over the RHS as per default
)

for pair in reverse(hds)
    h, delta = pair

    rainbow = RandbowDistribution(delta, h)

    disd = disorder_entropy_average(ps, n_trials, rainbow, 
        min_window_size = min_window_size,
        window_increment = window_increment,
        L_over_2 = L_over_2, 
        window_max = L_over_2,
        L = L,
        l_range = l_range,
        l_length = l_length,
        chain_start = chain_start)

    filter!(:l => x -> x  < maximum(disd.l) - 2, disd) 

    @df disd scatter!(fig_6, :l, :entropy, label = L"h = %$h, \delta = %$delta",
        xlabel = L"l", ylabel = L"S_A", clip_mode = "individual")

    # Fit square root 
    optim_result = optimize(x -> sqrt_loss(x, disd.l, disd.entropy), ones(2))
    a, b = Optim.minimizer(optim_result)
    plot!(fig_6, x -> a + (b*sqrt(x)), 10:0.1:100, 
    label = (pair == first(hds) ? "SDRG" : nothing),
    colour = :red, linestyle = :dashdot)

end

saveplot(fig_6, img_path, "randbow_ratio_entropy")