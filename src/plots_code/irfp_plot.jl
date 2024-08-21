# Build IRFP plot
println("Generating IFRP plot")
irfp_L = 5 * 10^4

irfp_params = SDRGParams(L = irfp_L, plots = true, boundaries = :periodic, 
    verbose = true, img_path = img_path, data_path = data_path, Δ = 1.0)

irfp_singlets = fill((0, 0), irfp_L ÷ 2)
irfp_Js = rand(uniform_distribution, irfp_L)
irfp_system = generate_system(irfp_Js, irfp_L)
sdrg!(irfp_singlets, irfp_system, irfp_params)
println("") #Not starting a new line properly?