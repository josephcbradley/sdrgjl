# P_r distribution
L = 1000
ps = SDRGParams(L = L,
    plots = false, 
    boundaries = :open, 
    verbose = false, 
    Δ = 0.0)
h = BigFloat(7.0)
randbow = RandbowDistribution(BigFloat(1.0), h)

bins = 2:2:100
hist = fit(Histogram, Float64[], bins, closed = :left)
@showprogress "Doing h = 7 P_r distribution..." for n in 1:10^4 
    singlets = fill((0, 0), L ÷ 2)
    Js = rand(randbow, L - 1)
    system = generate_system(Js, L)
    sdrg!(singlets, system, ps)
    coding = rainbow_coding(singlets)
    new_sequences = (count_consecutive_ones(coding))
    new_hist = fit(Histogram, new_sequences, bins)
    merge!(hist, new_hist)
end
normed_h_7 = normalize(hist, mode = :pdf)

h = BigFloat(10.0)
randbow = RandbowDistribution(BigFloat(1.0), h)

hist = fit(Histogram, Float64[], bins, closed = :left)
@showprogress "Doing h = 10 P_r distribution..." for n in 1:10^4 
    singlets = fill((0, 0), L ÷ 2)
    Js = rand(randbow, L - 1)
    system = generate_system(Js, L)
    sdrg!(singlets, system, ps)
    coding = rainbow_coding(singlets)
    new_sequences = (count_consecutive_ones(coding))
    new_hist = fit(Histogram, new_sequences, bins)
    merge!(hist, new_hist)
end

normed_h_10 = normalize(hist, mode = :pdf)

p_r_distribution_plot = plot(normed_h_10, dpi = 600, xlab = L"l_r",
    ylab = L"P_r", label = L"h = %$h", color = :orange,
    yscale = :log10, 
    ylims = (10^-6, 1),
    fillrange = 10^-6,
    legend = :topright
    )

plot!(normed_h_7, label = L"h = 7.0", color = :aqua, 
    yscale = :log10, 
    fillrange = 10^-6
    )

saveplot(p_r_distribution_plot, img_path, "p_r_distribution")