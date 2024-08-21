#Run everything!
println("Loading code...")
using Revise, Distributions, LaTeXStrings, StatsPlots, CSV, GLM, Optim, StatsBase
includet("../load.jl")

rainbow = RandbowDistribution{BigFloat}(0, 5)
L = 50
Js = rand(rainbow, L - 1)
system = generate_system(Js, L)
ps = SDRGParams(L, false, :open, false, "", "", 0.0)
singlets = fill((0, 0), L ÷ 2)
sdrg!(singlets, system, ps)
min_window_size = 1
window_increment = 1
r = 10
start_pos = (L - r) ÷ 2 
window_adjust = -1
window_ceiling = (L - r) ÷ 2
buffer = generate_negativity_buffer(min_window_size, window_increment, window_ceiling)
negativity_spectrum!(
    buffer, singlets, r, 
    L = L, 
    min_window_size = min_window_size, 
    window_increment = window_increment, 
    window_adjust = window_adjust,
    start_pos = start_pos,
    window_ceiling = window_ceiling
)
buffer