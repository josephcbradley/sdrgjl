println("Running tests...")
using Test, Distributions, GLM

include("load.jl")
println("Code loaded")

@testset "SDRG Tests" begin
    # Regular SDRG test
    L = 4
    singlets = fill((0, 0), L ÷ 2)
    Js = Float64.([1, 4, 3, 2])
    system = generate_system(Js, L)
    ps = SDRGParams(L = L, plots = false, boundaries = :periodic, 
        verbose = false, Δ = 1.0)
    sdrg!(singlets, system, ps)
    @test singlets == [(2, 3), (1, 4)]

    #Right hand boundary test 
    L = 6
    singlets = fill((0, 0), L ÷ 2)
    Js = Float64.([1, 4, 3, 2, 5])
    system = generate_system(Js, L)
    ps = SDRGParams(L = L, plots = false, 
        boundaries = :open, 
        verbose = false, Δ = 1.0)
    sdrg!(singlets, system, ps)
    @test singlets == [(5, 6), (2, 3), (1, 4)]

    #Left hand boundary test 
    L = 6
    singlets = fill((0, 0), L ÷ 2)
    Js = Float64.([5, 2, 3, 4, 1])
    system = generate_system(Js, L)
    ps = SDRGParams(L = L, plots = false, 
        boundaries = :open, 
        verbose = false, Δ = 1.0)
    sdrg!(singlets, system, ps)
    @test singlets == [(1, 2), (4, 5), (3, 6)]

    #rainbow test
    L = 6
    singlets = fill((0, 0), L ÷ 2)
    rainbow = RandbowDistribution(0.0, 100.0)
    Js = rand(rainbow, L - 1)
    system = generate_system(Js, L)
    ps = SDRGParams(L = L, plots = false, boundaries = :open, verbose = false,
    Δ = 1.0)
    sdrg!(singlets, system, ps)
    @test singlets == [(3, 4), (2, 5), (1, 6)]
end
    

@testset "Region Tests" begin
    @test region_test(52, 51, 60) == true
    @test region_test(50, 51, 60) == false
    @test region_test(51, 51, 60) == true
end

@testset "Counting tests" begin 
    singlets = [(1, 2), (3, 4), (5, 6), (7, 8)]
    @test complement_A_count(singlets, (1, 4)) == 0
    @test complement_A_count(singlets, (4, 5)) == 2
    @test two_intervals_count(singlets, 1, 2, 0) == 0
    @test two_intervals_count(singlets, 1, 3, 0) == 1
    @test two_intervals_count(singlets, 1, 2, 1) == 0
    # Try getting something with an r != 1
    singlets = [(1, 2), (3, 5), (4, 6), (7, 8)]
    @test two_intervals_count(singlets, 1, 3, 1) == 1
end

@testset "Entropy tests" begin 
    #Simple four site entropy test 
    singlets = [(1, 2), (3, 4)]
    A = (1, 2)
    @test entanglement_entropy(singlets, A) == 0.0
    A = (2, 3)
    @test entanglement_entropy(singlets, A) == 2 * log(2)
    #Longer entropy test
    singlets = ([(i % 20 + 1, (i + 1) % 20 + 1) for i in 0:20-1])
    A = (11, 20)
    @test entanglement_entropy(singlets, A) == 2 * log(2)
    
    # Calculate coefficient of entropy plot
    distribution = Uniform(0, 1)
    params_1000 = SDRGParams(L = 1000, plots = false, boundaries = :periodic,
        Δ = 1.0)
    println("Testing entropy disorder average.")
    results_1000 = disorder_entropy_average(params_1000, 2_000, distribution)
    adjust_l_entropy_data!(results_1000, params_1000.L)
    results_1000.l_adj_log10 = log.(Ref(10), results_1000.l_adj)
    ols = lm(@formula(entropy ~ l_adj_log10), results_1000)
    k = coef(ols)[2]
    @test isapprox(k, 0.5, atol = 0.1)
end

@testset "Negativity tests" begin 
    singlets = [(1, 15), (2, 19), (3, 4), (12, 16)]
    A₁ = (1, 10)
    l = 10
    r = 0
    @test negativity(singlets, A₁[1], l, r) == 2 * log(2)
    L = 1000
    ps = SDRGParams(L = L, boundaries = :open)
    rainbow = RandbowDistribution{BigFloat}(0, 5)
    Js = rand(rainbow, L - 1)
    system = generate_system(Js, L)
    singlets = fill((0, 0), L ÷ 2)
    sdrg!(singlets, system, ps)
    start_pos = L ÷ 2
    min_window_size = 1
    window_increment = 1
    window_ceiling = L ÷ 2
    window_adjust = -1
    buffer = generate_negativity_buffer(
        min_window_size,
        window_increment,
        window_ceiling
    )
    negativity_spectrum!(buffer, singlets, 0,
        L = L, min_window_size = min_window_size, 
        window_increment = window_increment, 
        window_ceiling = window_ceiling, 
        start_pos = start_pos, 
        window_adjust = window_adjust)
    @test buffer.negativity ./ log(2) ≈ buffer.l
end
