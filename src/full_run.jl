#Run everything!
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
n_trials = 50_000
println("Using $n_trials trials for all disorder averages.")

pgfplotsx()
include(plots_code_path * "irfp_plot.jl")
include(plots_code_path * "disorder_entropy_baseline_pbc.jl")
include(plots_code_path * "disorder_entropy_baseline_obc.jl")
include(plots_code_path * "disorder_negativity_baseline_pbc.jl")
include(plots_code_path * "disorder_negativity_baseline_obc.jl")
include(plots_code_path * "randbow_sqrt_entropy.jl")
include(plots_code_path * "randbow_ratio_entropy.jl")
include(plots_code_path * "pr_distribution.jl")
include(plots_code_path * "pb_distribution.jl")
include(plots_code_path * "randbow_central_negativity_obc_XX.jl")
include(plots_code_path * "randbow_central_negativity_obc_XXY.jl")
include(plots_code_path * "randbow_pure_negativity.jl")
include(plots_code_path * "powerbow_entropy.jl")
include(plots_code_path * "powerbow_central_negativity_obc_XX.jl")
include(plots_code_path * "powerbow_central_negativity_obc_XXY.jl")
include(plots_code_path * "exact_solve_randbow.jl")
include(plots_code_path * "exact_solve_randbow_clean.jl")
include(plots_code_path * "exact_solve_powerbow.jl")
include(plots_code_path * "exact_solve_powerbow_clean.jl")
include(plots_code_path * "randbow_central_neg_XX_vary_r.jl")
include(plots_code_path * "randbow_central_neg_XXY_vary_r.jl")
include(plots_code_path * "powerbow_central_neg_XX_vary_r.jl")
include(plots_code_path * "powerbow_central_neg_XXY_vary_r.jl")
include(plots_code_path * "timings_sdrg.jl")
include(plots_code_path * "timings_analysis.jl")




