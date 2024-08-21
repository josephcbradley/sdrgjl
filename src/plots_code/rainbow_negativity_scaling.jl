rainbow = RandbowDistribution{BigFloat}(0, 5)
L = 1000
Js = rand(rainbow, L - 1)
system = generate_system(Js, L)
ps = SDRGParams(L, false, :open, false, "", "", 0.0)
singlets = fill((0, 0), L ÷ 2)
sdrg!(singlets, system, ps)
min_window_size = 1
window_increment = 1
start_pos = L ÷ 2 
window_adjust = -1
r = 0
window_ceiling = (L - r) ÷ 2
buffer = generate_negativity_buffer(min_window_size, window_increment, window_ceiling)
negativity_spectrum!(
    buffer, singlets, 0, 
    L = L, 
    min_window_size = min_window_size, 
    window_increment = window_increment, 
    window_adjust = window_adjust,
    start_pos = start_pos,
    window_ceiling = window_ceiling
)
buffer.count = buffer.negativity ./ log(2)
@df buffer scatter(:l, :negativity)

saveplot(p_r_distribution_plot, img_path, "rainbow_negativity_scaling")