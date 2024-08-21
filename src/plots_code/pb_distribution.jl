# Now bubble dsitribution 
L = 1000
bins = 1:2:100
h = BigFloat(7.0)
randbow = RandbowDistribution(BigFloat(1.0), h)

hist = fit(Histogram, Float64[], bins, closed = :left)
@showprogress "Doing h = 7 P_b distribution..." for n in 1:10^4
    singlets = fill((0, 0), L ÷ 2)
    Js = rand(randbow, L - 1)
    system = generate_system(Js, L)
    sdrg!(singlets, system, ps)
    coding = rainbow_coding(singlets)
    new_sequences = (count_consecutive_ones(coding, sign = false))
    new_hist = fit(Histogram, new_sequences, bins)
    merge!(hist, new_hist)
end
normed_h_7 = normalize(hist, mode = :pdf)

h = BigFloat(10.0)
randbow = RandbowDistribution(BigFloat(1.0), h)

hist = fit(Histogram, Float64[], bins, closed = :left)
@showprogress "Doing h = 10 P_b distribution..." for n in 1:10^4 
    singlets = fill((0, 0), L ÷ 2)
    Js = rand(randbow, L - 1)
    system = generate_system(Js, L)
    sdrg!(singlets, system, ps)
    coding = rainbow_coding(singlets)
    new_sequences = (count_consecutive_ones(coding, sign = false))
    new_hist = fit(Histogram, new_sequences, bins)
    merge!(hist, new_hist)
end
normed_h_10 = normalize(hist, mode = :pdf)

p_b_distribution_plot = plot(normed_h_10, dpi = 600, xlab = L"l_b",
    ylab = L"P_b", label = L"h = %$h", color = :orange,
    yscale = :log10, xscale = :log10, xlims = (1, 100), 
    ylims = (10^-4, 1), fillrange = 10^-4,
    legend = :topright)

plot!(normed_h_7, label = L"h = 7.0", color = :aqua, yscale = :log10, xscale = :log10,
    fillrange = 10^-4)

saveplot(p_b_distribution_plot, img_path, "p_b_distribution")
