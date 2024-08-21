L = 2000
L_over_2 = L ÷ 2
ps = SDRGParams(L, false, :open, false, "", "", big(0.0))
αs = (BigFloat.([10, 100, 500]))
chain_start = L_over_2 + 1
min_window_size = 2 
window_increment = 10
l_range = min_window_size:window_increment:L_over_2
l_length = length(l_range)

fig_4 = plot(dpi = 600, xlabel = L"l", ylabel = L"S_A",
    #xlims = (2, 100 + 10), 
    #xscale = :log2, 
    #yscale = :log2, 
    legend = :bottomright)

for α in reverse(αs)
#α = 0.1

    powerbow = PowerbowDistribution(BigFloat(1.0), BigFloat(α))

    disd = disorder_entropy_average(ps, n_trials, powerbow, 
        min_window_size = min_window_size,
        window_increment = window_increment,
        L_over_2 = L_over_2, 
        window_max = L_over_2,
        L = L,
        l_range = l_range,
        l_length = l_length,
        chain_start = chain_start)
 

    @df disd scatter!(fig_4, :l, :entropy, label = L"\alpha = %$α")

    #= # Fit square root 
    # Not coverging for low L, so only fit from 
    # l = 50:l_max - 10
    disd_filter = @view disd[30:end - 10, :]
    optim_result = optimize(x -> sqrt_loss(x, disd_filter.l, disd_filter.entropy), ones(2))
    a, b = Optim.minimizer(optim_result)
    plot!(fig_4, x -> a + (b*sqrt(x)), 2:0.1:L_over_2 - 10, label = nothing,
    colour = :red, linestyle = :dashdot) =#

end
fig_4

saveplot(fig_4, img_path, "powerbow_entropy_strong_h_2000")

#= L = 1000
L_over_2 = L ÷ 2
ps = SDRGParams(L, false, :open, false, "", "", 0.0)
chain_start = L_over_2 + 1
min_window_size = 2 
window_increment = 10
l_range = min_window_size:window_increment:L_over_2
l_length = length(l_range)

fig_4 = plot(dpi = 600, xlabel = L"l", ylabel = L"S_A",
    #xlims = (2, 100 + 10), 
    #xscale = :log2, 
    #yscale = :log2, 
    legend = :bottomright)

for α in reverse(αs)
#α = 0.1

    powerbow = PowerbowDistribution(BigFloat(1.0), BigFloat(α))

    disd = disorder_entropy_average(ps, n_trials, powerbow, 
        min_window_size = min_window_size,
        window_increment = window_increment,
        L_over_2 = L_over_2, 
        window_max = L_over_2,
        L = L,
        l_range = l_range,
        l_length = l_length,
        chain_start = chain_start)
 

    @df disd scatter!(fig_4, :l, :entropy, label = L"α = %$α")

    #= # Fit square root 
    # Not coverging for low L, so only fit from 
    # l = 50:l_max - 10
    disd_filter = @view disd[30:end - 10, :]
    optim_result = optimize(x -> sqrt_loss(x, disd_filter.l, disd_filter.entropy), ones(2))
    a, b = Optim.minimizer(optim_result)
    plot!(fig_4, x -> a + (b*sqrt(x)), 2:0.1:L_over_2 - 10, label = nothing,
    colour = :red, linestyle = :dashdot) =#

end
fig_4

saveplot(fig_4, img_path, "powerbow_entropy_strong_h_1000") =#