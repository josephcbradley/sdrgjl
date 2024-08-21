using BenchmarkTools, Plots, StatsPlots, BenchmarkPlots, Revise
includet("load.jl")

# Region benchmarks 
region_test_suite = BenchmarkGroup()
region_test_suite["left_region_test"] = @benchmarkable left_region_test($52, $51, $60, $100)
region_test_suite["new_left_region_test"] = @benchmarkable better_left($52, $51, $60)
region_test_suite["right_region_test"] = @benchmarkable right_region_test($26, $10, $21, $30, $100)
tune!(region_test_suite)
region_results = run(region_test_suite)
plot(region_results, ylims = (0, 10))

# Analysis benchmarks 
analysis_test_suite = BenchmarkGroup()
analysis_test_suite["entropy"] = @benchmarkable singlet_entropy($([(1, 15), (2, 19), (3, 4), (12, 16)]), $((1, 50)), $100)
analysis_test_suite["new_entropy"] = @benchmarkable singlet_entropy_new($([(1, 15), (2, 19), (3, 4), (12, 16)]), $((1, 50)), $100)
tune!(analysis_test_suite)
analysis_results = run(analysis_test_suite)
plot(analysis_results)

# Entropy benchmarks
entropy_test_suite = BenchmarkGroup()
#Pre
using Distributions
uniform_distribution = Uniform(0, 1)
L = 2000
ps = SDRGParams(L, false)
L_over_2 = L ÷ 2
links = WrappedVec([(i % L + 1, (i + 1) % L + 1) for i in 0:L-1])
singlets = generate_empty_singlets()
Js = rand(uniform_distribution, L)
system = System(deepcopy(links), WrappedVec(Js))
buffer = generate_buffer()
sdrg!(singlets, system, buffer, ps)
entropy_test_suite["old"] = @benchmarkable entropy_spectrum!(df, $singlets.links) setup =(df = DataFrame(l = collect(10:10:L_over_2), entropy = zeros(length(10:10:L_over_2))))
#entropy_test_suite["new"] = @benchmarkable entropy_spectrum_new!(df, $singlets.links) setup =(df = DataFrame(l = collect(10:10:L_over_2), entropy = zeros(length(10:10:L_over_2))))
entropy_test_suite["new_sliding"] = @benchmarkable entropy_spectrum_new_sliding!(df, $singlets.links) setup =(df = DataFrame(l = collect(10:10:L_over_2), entropy = zeros(length(10:10:L_over_2))))
tune!(entropy_test_suite)
ent_results = run(entropy_test_suite)
plot(ent_results)

buf = DataFrame(l = collect(10:10:L_over_2), entropy = zeros(length(10:10:L_over_2)))
@profview entropy_spectrum_new!(buf, singlets.links)