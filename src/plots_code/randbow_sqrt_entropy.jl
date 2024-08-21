println("Randbow sqrt entropy plot")
L = 240
L_over_2 = L ÷ 2
ps = SDRGParams(L, false, :open, false, "", "", 0.0)
hs = (BigFloat.([0.5, 3., 7., 10.]))
chain_start = L_over_2 + 1
min_window_size = 2
window_increment = 1
l_range = min_window_size:window_increment:L_over_2
l_length = length(l_range)

fig_4 = plot(dpi = 600, xlabel = "l", ylabel = "S",
    xlims = (2, 100 + 10), 
    xscale = :log10, 
    yscale = :log10, legend = :topleft)

for h in reverse(hs)

    rainbow = RandbowDistribution(BigFloat(1.0), BigFloat(h))

    disd = disorder_entropy_average(ps, n_trials, rainbow, 
        min_window_size = min_window_size,
        window_increment = window_increment,
        L_over_2 = L_over_2, 
        window_max = L_over_2,
        L = L,
        l_range = l_range,
        l_length = l_length,
        chain_start = chain_start)

    #Fudge - drop the last couple of values of l 
    filter!(:l => x -> x  < maximum(disd.l) - 2, disd) 

    @df disd scatter!(fig_4, :l, :entropy, label = L"h = %$h", clip_mode = "individual")

    # Fit square root 
    optim_result = optimize(x -> sqrt_loss(x, disd.l, disd.entropy), ones(2))
    a, b = Optim.minimizer(optim_result)
    plot!(fig_4, x -> a + (b*sqrt(x)), 10:0.1:100, label = (h == first(hs) ? "SDRG" : nothing),
    colour = :red, linestyle = :dashdot)

end

saveplot(fig_4, img_path, "randbow_sqrt_entropy")