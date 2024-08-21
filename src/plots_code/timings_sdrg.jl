using DataFrames, Distributions, BenchmarkTools, LaTeXStrings
setprecision(BigFloat, 128)
data_types = [Float64, BigFloat]
Ls = [i for i in 500:500:10000]

data = DataFrame(data_type = String[], L = Int64[], median_time_mu_s = Float64[])
for pair in Iterators.product(data_types, Ls)
    data_type, L = pair
    singlets = fill((0, 0), L ÷ 2)
    ps = SDRGParams(L = L, boundaries = :periodic)
    uniform = Uniform{data_type}(0, 1)
    Js = rand(uniform, L)
    system = generate_system(Js, L)
    res = @benchmark sdrg!($singlets, $system, $ps) evals = 10
    push!(data, (string(data_type), L, median(res).time / 1000.))
end

replace!(data.data_type, "BigFloat" => "128-bit", "Float64" => "64-bit")

plt = @df data plot(:L, :median_time_mu_s, group = :data_type,
    xlabel = L"L", ylabel = L"\mu\textrm{s}", xlims = (0, 10000), ylims = (3, 11));

saveplot(plt, img_path, "sdrg_timings");
