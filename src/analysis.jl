using Accessors, DataFrames, ProgressMeter, OnlineStats, Polyester

export entanglement_entropy, entropy_spectrum!, Y, Lc, adjust_l_entropy_data!, n_means, disorder_entropy_average, shift_negativity, negativity, negativity_spectrum!, disorder_negativity_average, add_shifted_negativity!, disorder_negativity_average, generate_analysis_buffer, generate_entropy_buffer, generate_negativity_buffer, singlets_to_row_form, rainbow_region_counting

function entanglement_entropy(singlets, A)::Float64
    # The entanglement entropy Sₐ is obtained by determining 
    # the number of singlets between A and B for each realisation of 
    # the disorder and then by multiplying the average by ln(2). 

    # Calculate the entropy for one window position of one realisation
    
    complement_A_count(singlets, A) * log(2)

end

function entanglement_entropy!(vec, pos, singlets, A)::Float64
    # The entanglement entropy Sₐ is obtained by determining 
    # the number of singlets between A and B for each realisation of 
    # the disorder and then by multiplying the average by ln(2). 

    # Calculate the entropy for one window position of one realisation
    
    vec[pos] = complement_A_count(singlets, A) * log(2)

end

#= function entanglement_entropy(singlets, window::R) where R <: AbstractRange{Int64}
    singlet_entropy(singlets, (window.start, window.stop))
end =#

function entropy_spectrum!(buffer, singlets;
    min_window_size = 10,
    window_increment = 10,
    L_over_2 = length(singlets), 
    L = L_over_2 * 2,
    l_range = min_window_size:window_increment:L_over_2,
    l_length = length(l_range),
    chain_start = 1)

    # the buffer should have columns l, and entropy
    if names(buffer) != ["l", "entropy"]
        error("Column names not as expected in entropy buffer!")
    end

    #the buffer should have l spaced from 10 to L / 2
    if size(buffer, 1) != l_length
        error("Buffer not long enough! Buffer size is $(size(buffer, 1)), l is $l_length")
    end

    # To calculate, we will take a window of length l, and then 
    # record it's entropy as we move the window along the singlets
    # Once we are done with a length l, we will increment the window 
    # size until it is at its maximum (L_over_2) and then stop
    # Start from the chain centre
    A = (chain_start, 1)
    # Just expand window from the left
    for i in eachindex(l_range)
        l = l_range[i]
        A = @set A[2] = A[1] + l - 1
        # To do averages, we note that for each window size l 
        # there are L - l possible window positions
        entanglement_entropy!(buffer.entropy, i, singlets, A)
    end

    return nothing
end

function Y(x)
    sin(x) * (1.0 + (4/3 * 0.115 * sin(x)^2))
end

function Lc(l, L)
    (L / π) * Y((π*l) / L)
end

function adjust_l_entropy_data!(df, L)
    #Check that columns we need are in the df
    if "l" ∉ names(df) || "entropy" ∉ names(df)
        error("Column names not as expected in entropy data!")
    end

    df.l_adj .= Lc.(df.l, Ref(L))
end

#= Disorder averaged entropy 
To calculate the disorder average of many instances of the system, 
we will set up a 'master dataframe' that will gather the results of 
disorder realsiations
=#

function n_means(N)
    # Create a vector of N means that can be fitted OnlineStats
    # and thread sage
    start = "Group("
    mean = "Mean(), "
    stop = ")"
    eval(Meta.parse(start * mean^N * stop))
end

function disorder_entropy_average(params, ntrials, bond_distribution;
    min_window_size = 10,
    window_increment = 10,
    L_over_2 = params.L ÷ 2, 
    window_max = L_over_2,
    L = params.L,
    l_range = min_window_size:window_increment:window_max,
    l_length = length(l_range),
    chain_start = 1)
    
    bond_type = rand(bond_distribution, 1) |> eltype

    #Create entropy vector 
    S = n_means(l_length)
    # Start progress meter
    p = Progress(ntrials, "Disordered entropy progress:")
    
    boundaries = params.boundaries

    Threads.@threads for trial in 1:ntrials
        # Singlets output buffer
        singlets = fill((0, 0), L ÷ 2)

        # bonds buffer - adjust for boundary conditions
        # tricky - sdrg procedure will push a zero if we have 
        # an odd number of Js, and then broadcasting will fail, 
        # so preempt that here by just making it equal to L
        Js = Vector{bond_type}(undef, L)

        #entropy buffer - also to rewrite every time
        entropy_buffer = DataFrame(l = collect(l_range), entropy = zeros(l_length))

        if boundaries == :open
            # carefully w/ the boundary conditions
            Js[1:end - 1] .= rand(bond_distribution, L - 1)
            Js[end] = zero(bond_type)
        elseif boundaries == :periodic
            Js .= rand(bond_distribution, L)
        end

        system = generate_system(Js, L)

        
        sdrg!(singlets, system, params)
        entropy_spectrum!(entropy_buffer, singlets,
        min_window_size = min_window_size,
        window_increment = window_increment,
        L_over_2 = L_over_2, 
        L = L, chain_start = chain_start,
        l_range = l_range,
        l_length = l_length)
        #merge results with master 
        fit!(S, entropy_buffer.entropy)
        next!(p)
    end

    #Allocate output
    out = DataFrame(l = collect(l_range), entropy = Vector{Float64}(undef, l_length))
    #Parse final values from s
    for i in eachindex(l_range)
        out.entropy[i] = S[i].μ
    end

    return out
end


function negativity(singlets, i, l, r)::Float64
    two_intervals_count(singlets, i, l, r) * log(2)
end

function negativity!(negativity_vec, pos, singlets, i, l, r)::Float64
    negativity_vec[pos] = two_intervals_count(singlets, i, l, r) * log(2)
end


function negativity_spectrum!(buffer, singlets, r;
    L = length(singlets)  * 2,
    min_window_size = 10,
    window_increment = 10,
    window_ceiling = round(Int64, (L - r)/2, RoundDown),
    window_adjust = 0,
    l_range = min_window_size:window_increment:window_ceiling,
    l_length = length(l_range),
    start_pos = 1)

    # the buffer should have columns l, and negativity
    if names(buffer) != ["l", "negativity"]
        error("Column names not as expected in negativity buffer!")
    end

    #the buffer should have l spaced from 10 to L / 2
    if size(buffer, 1) != l_length
        error("Buffer not long enough! Buffer size is $(size(buffer, 1)), l_length is $l_length")
    end

    # To calculate, we will take a window of length l, and then 
    # record its negativity as we expand the window along the singlets
    # Once we are done with a length l, we will increment the window 
    # size until it is at its maximum (L_over_2) and then stop

    #Record position that we can adjust with window_adjust
    position = start_pos

    for i in 1:l_length
        l = l_range[i]
        position = start_pos + ((i - 1) * window_adjust)
        negativity!(buffer.negativity, i, singlets, position, l, r)
    end

    return nothing

end

function disorder_negativity_average(params, ntrials, bond_distribution, r;
    start_pos = 1,
    min_window_size = 10,
    window_increment = 10,
    window_adjust = 0,
    L = params.L,
    window_ceiling = round(Int64, (L - r)/2, RoundDown),
    l_range = min_window_size:window_increment:window_ceiling,
    l_length = length(l_range))

    #Create negativity vector 
    ϵ = n_means(l_length)
    #Links will be the same every time, so no need to re calculate them
    p = Progress(ntrials, "Disordered negativity progress:")
    # Singlets output buffer
    boundaries = params.boundaries
    # Determine bond type 
    bond_type = eltype(rand(bond_distribution, 1))

    Threads.@threads for trial in 1:ntrials
        # bonds buffer - adjust for boundary conditions
        # tricky - sdrg procedure will push a zero if we have 
        # an odd number of Js, and then broadcasting will fail, 
        # so preempt that here by just making it equal to L
        Js = Vector{bond_type}(undef, L)
        if boundaries == :open
            # carefully w/ the boundary conditions
            Js[1:end - 1] .= rand(bond_distribution, L - 1)
            Js[end] = zero(bond_type)
        elseif boundaries == :periodic
            Js .= rand(bond_distribution, L)
        end

        singlets = fill((0, 0), L ÷ 2)
        system = generate_system(Js, L)
        sdrg!(singlets, system, params)
        negativity_buffer = DataFrame(l = collect(l_range), negativity = zeros(l_length))
        negativity_spectrum!(negativity_buffer, singlets, r,
        L = L, min_window_size = min_window_size, 
        window_increment = window_increment,
        window_adjust = window_adjust,
        window_ceiling = window_ceiling, l_range = l_range,
        l_length = l_length, start_pos = start_pos)
        #merge results with master 
        fit!(ϵ, negativity_buffer.negativity)
        next!(p)
    end

    #Allocate output
    out = DataFrame(l = collect(l_range), negativity = Vector{Float64}(undef, l_length))
    #Parse final values from s
    for i in eachindex(l_range)
        out.negativity[i] = ϵ[i].μ
    end

    return out
end

function shift_negativity(ϵ, L)
    ϵ - ((log(2) / 6)*log(L))
end

function add_shifted_negativity!(df, L)
    #Check that columns we need are in the df
    if "l" ∉ names(df) || "negativity" ∉ names(df)
        error("Column names not as expected in negativity data!")
    end

    df.shifted_negativity .= shift_negativity.(df.negativity, Ref(L))

end

function generate_analysis_buffer(variable,
    min_window_size,
    window_increment,
    window_ceiling;
    l_range = min_window_size:window_increment:window_ceiling,
    l_length = length(l_range))

    DataFrame("l" => collect(l_range), String(variable) => Vector{Float64}(undef, l_length))
end

function generate_entropy_buffer(
    min_window_size,
    window_increment,
    window_ceiling;
    l_range = min_window_size:window_increment:window_ceiling,
    l_length = length(l_range))

    generate_analysis_buffer(:entropy, 
    min_window_size,
    window_increment,
    window_ceiling;
    l_range = min_window_size:window_increment:window_ceiling,
    l_length = length(l_range))
end

function generate_negativity_buffer(
    min_window_size,
    window_increment,
    window_ceiling;
    l_range = min_window_size:window_increment:window_ceiling,
    l_length = length(l_range))

    generate_analysis_buffer(:negativity, 
    min_window_size,
    window_increment,
    window_ceiling;
    l_range = min_window_size:window_increment:window_ceiling,
    l_length = length(l_range))
end


#= ECDF of rainbow sites 

Given a set of singlets S, how to determine the 
distribution function of rainbow sections? 

Start by creating a function to conert individual 
    singlet data into the 'row representation' a la 
    Paola.
=#

function singlets_to_row_form(singlets)
    N = length(singlets)
    L = 2N
    out = Vector{Int64}(undef, L)
    for i in 1:N
        @inbounds l, r = singlets[i]
        @inbounds out[l] = i
        @inbounds out[r] = i
    end
    out
end

#= Want to estimate a PDF of rainbow region lengths
Start by counting how many regions of different lengths there are
=#

function rainbow_region_counting(singlet_row::Vector{Int64})
    #assume we have singlets in singlet row form 
    L = length(singlet_row)
    #longest possible rainbow section would be L ÷ 2
    data = Vector{Int64}()
    current_streak = 0
    i = 1
    while i <= L ÷ 2
        # if rainbow streak continues 
        current_site = singlet_row[i]
        next_site = singlet_row[i + 1]
        if current_site != next_site
            current_streak += 1
            i += 1
        else
            push!(data, current_streak)
            current_streak = 0
            i += 2
        end
    end
    filter!(x -> x != 0, data)
    data .+ 1
end

function rainbow_coding(singlets::Vector{Tuple{Int64, Int64}})
    L_over_2 = length(singlets)
    L = L_over_2 * 2
    #sort!(singlets, by = x -> x[1])
    coding = Vector{Bool}(undef, L)
    #code 0 for bubble, 1 for rainbow

    for singlet in singlets 
        a, b = singlet 
        if abs(a - b) == 1 #if they are only one apart, it's a bubble
            coding[a] = false
            coding[b] = false
        else
            coding[a] = true
            coding[b] = true
        end 
    end

    #assume that bridge link is a rainbow
    coding[L_over_2] = true
    #return the half chain 
    coding[1:L_over_2]
end

function region_counts(coding::Vector{Bool})
    #Assume that pmf is too long 
    L = length(coding)
    half = @view coding[begin:L ÷ 2]    
    counts = Int64[]
    curr_count = 0
    i = 1
    while i <= L ÷ 2
        if half[i] == 1
            peek = i + 1
            count = 1
            while half[peek]
                peek += 1
                count += 1
            end
            push!(counts, count)
            i = peek + 1
        end
        i + 1
    end
    count
end


#For fitting sqrt functions 
function sqrt_loss(x, input_range, data)
    a, b = x
    predict = a .+ (b .* sqrt.(input_range))
    sum((data .- predict) .^ 2)
end

function negativity_loss(x, input_range, data, L)
    k = x[1]
    predict = negativity_exact_pbc.(input_range, Ref(L), Ref(k))
    sum((data .- predict) .^ 2)
end

function count_consecutive_ones(coding; sign = true)
    out = Int64[]
    count = 0
    for x in coding 
        if x == sign
            count += 1
        else 
            push!(out, count)
            count = 0 
        end 
    end 

    #update with last count if we ended ona  rainbow
    push!(out, count)
    filter!(x -> x != zero(eltype(x)), out)
    out
    
end

function negativity_exact_pbc(l, L, k)
    #Note that if you want to plot this on a shifted neg graph, 
    # you need L = 1 and pass in l ∈ 0.0:0.5 
    top = Y(l*π / L)^2
    bottom = Y(2l*π / L)
    return (log(2)/6 * log(top / bottom)) + k
end

function generate_negativity_r_buffer(r_min, r_step, r_max)
    r_range = r_min:r_step:r_max
    return DataFrame(
        r = collect(r_range),
        negativity = Vector{Float64}(undef, length(r_range))
    )
end

function negativity_spectrum_over_r!(buffer, l, singlets;
    L = 1000)

    r_step = buffer.r[2] - buffer.r[1]

    # Check that r is a multiple of two 
    if r_step % 2 != 0
        error("r step is not a multiple of two! It is $(r_step).")
    end

    # the buffer should have columns r, and negativity
    if names(buffer) != ["r", "negativity"]
        error("Column names not as expected in negativity buffer!")
    end

    window_adjust = -(r_step ÷ 2)
    start_pos = (L ÷ 2) + window_adjust - l + 1

    #l is fixed and r varies

    for i in eachindex(buffer.r)
        r = buffer.r[i]
        #NB position is the LH edge of the LH interval 
        # window_adjust is just the r_step over two 
        position = start_pos + ((i - 1) * window_adjust)
        negativity!(buffer.negativity, i, singlets, position, l, r)
        #negativity(singlets, position, l, r) / log(2)
    end

    return nothing
end

function disorder_negativity_r_average(params, n_trials, bond_distribution, l;
    min_window_size = 2,
    window_increment = 2,
    L = params.L,
    window_ceiling = 10,
    r_range = min_window_size:window_increment:window_ceiling,
    r_length = length(r_range))

    #Create negativity vector 
    ϵ = n_means(r_length)
    #Links will be the same every time, so no need to re calculate them
    p = Progress(n_trials, "Disordered negativity r progress:")
    # Singlets output buffer
    boundaries = params.boundaries
    # Determine bond type 
    bond_type = eltype(rand(bond_distribution, 1))

    Threads.@threads for trial in 1:n_trials
        # bonds buffer - adjust for boundary conditions
        # tricky - sdrg procedure will push a zero if we have 
        # an odd number of Js, and then broadcasting will fail, 
        # so preempt that here by just making it equal to L
        Js = Vector{bond_type}(undef, L)
        if boundaries == :open
            # carefully w/ the boundary conditions
            Js[1:end - 1] .= rand(bond_distribution, L - 1)
            Js[end] = zero(bond_type)
        elseif boundaries == :periodic
            Js .= rand(bond_distribution, L)
        end

        singlets = fill((0, 0), L ÷ 2)
        system = generate_system(Js, L)
        sdrg!(singlets, system, params)
        negativity_buffer = DataFrame(r = collect(r_range), negativity = zeros(r_length))
        negativity_spectrum_over_r!(negativity_buffer, l, singlets, L = L)
        #merge results with master 
        fit!(ϵ, negativity_buffer.negativity)
        next!(p)
    end

    #Allocate output
    out = DataFrame(l = collect(r_range), negativity = Vector{Float64}(undef, r_length))
    #Parse final values from s
    for i in eachindex(r_range)
        out.negativity[i] = ϵ[i].μ
    end

    return out
end

function infinite_chain_negativity(l, k)
    (log(2) / 6 * log(
        l^2 / 2l
    )) + k
end

function infinite_chain_negativity_loss(x, input_range, data)
    k = x[1]
    predict = infinite_chain_negativity.(input_range, Ref(k))
    sum((data .- predict) .^ 2)
end

function cij(i, j, vectors)
    out = zero(eltype(vectors))
    for ϕ in eachcol(vectors)
        out += ϕ[i]'*ϕ[j]
    end
    out
end

function eig_entropy(λs)
    out = zero(eltype(λs))
    for λ in λs
        if λ < 1 && λ > 0
            out -= λ * log(λ) + (1 - λ)*log(1 - λ)
        end
    end
    out
end