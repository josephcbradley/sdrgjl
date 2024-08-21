#Run only some code 
println("Loading code...")
using Revise, Distributions, LaTeXStrings, Plots, StatsPlots, CSV, GLM, Optim, StatsBase
const font_size = 18
default(tickfontsize = font_size, yguidefontsize = font_size, xguidefontsize = font_size, legendfontsize = font_size, 
    widen = true)
includet("load.jl")

# Run tests
import InteractiveUtils: versioninfo
println("There are $(length(Base.loaded_modules)) loaded modules.")
versioninfo()

const img_path = "data/imgs/"
const data_path = "data/"
const plots_code_path = "plots_code/"
const file_types = [".png", ".tikz"]
const uniform_distribution = Uniform{BigFloat}(0, 1)
n_trials = 10_000
println("Using $n_trials trials for all disorder averages.")

pgfplotsx()
include(plots_code_path * "powerbow_entropy.jl")




