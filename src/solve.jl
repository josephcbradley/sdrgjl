#= System solvers 
=#

export sdrg!, J′_calc, Γ_calc, IRFP_density, β_calc, chain_wrap

using LinearAlgebra, Dictionaries, Plots, Accessors, ProgressMeter, UnPack, StatsBase, LaTeXStrings
pgfplotsx()

function J′_calc(Jₗ, Jᵣ, Jₘ, Δ)
    (Jₗ * Jᵣ) / ( (1.0 + Δ) * Jₘ)
end



function Γ_calc(J₀, Jₘ)::Real
    log(J₀ / Jₘ)
end

function IRFP_density(β, Γ)::Real
    (1 / Γ) * exp(-β / Γ)
end

function β_calc(Jₘ, Js)::Vector{Real}
    log.(Ref(Jₘ) ./ Js)
end

function chain_wrap(x, L)
    if x < oneunit(x)
        L - x
    elseif x > L
        x - L
    else
        x
    end
end

function sdrg!(singlets, system, params)
    pgfplotsx()
    @debug "Bond type is $(eltype(system.bonds))"
    if params.boundaries == :periodic
        sdrg_periodic!(singlets, system, params)
    elseif params.boundaries == :open
        sdrg_open!(singlets, system, params)
    else
        error("Boundary type not recognised!")
    end
end

function sdrg_periodic!(singlets, system, params)
    #Overwite singlets with output of SDRG process

    # Unpack parameters and data
    @unpack L, plots, boundaries, verbose, img_path, data_path, Δ = params
    @unpack nodes, active, bonds = system

    #Initialise plot
    if plots
        plt = plot(dpi = 300, ylims = (0, 1), xlims = (0, 15), legend= :bottomright)
        hist_points = (round.(Ref(Int64), L/2 .* [0.1, 0.5, 0.7, 0.99]))
    end

    #main loop
    if verbose 
        #println("Bond type is $(eltype(bonds)).")
        p = Progress(L ÷ 2, "IRFP Progress:")
    end

    n_active = sum(active)
    m = 1

    # Make a note of J₀ if we need it for plots

    if plots 
        J₀ = maximum(bonds)
    end


    while n_active > 2
        curr_bonds = @view bonds[active]
        curr_L = length(curr_bonds) # should this be n_active? don't think so
        curr_nodes = @view nodes[active]

        Jₘ, Jₘ_position = findmax(curr_bonds)

        # Calculate adjacent spins
        Jₗ_position = chain_wrap(Jₘ_position - 1, curr_L)
        Jᵣ_position = chain_wrap(Jₘ_position + 1, curr_L)
        Jₗ = curr_bonds[Jₗ_position]
        Jᵣ = curr_bonds[Jᵣ_position]

        J_prime = J′_calc(Jₗ, Jᵣ, Jₘ, Δ)

        # Set Jₗ to J_prime
        curr_bonds[Jₗ_position] = J_prime

        #Deactivate sites at Jₘ, Jᵣ
        curr_active = @view active[active]
        curr_active[Jₘ_position] = false
        curr_active[Jᵣ_position] = false

        #add singlet 
        singlets[m] = (curr_nodes[Jₘ_position], curr_nodes[Jᵣ_position])

        m += 1

        n_active = sum(active)

        if plots && in(m, hist_points)
            β = β_calc(maximum(curr_bonds), curr_bonds)
            histogram!(plt, view(β, β .<= 15.0), label = L"m = %$m", 
            normalize = true, bins = 30, xlims = (0, 15))


            if m == last(hist_points)
                #println("Adding IRFP inset...")
                histogram!(plt,
                    β, label = L"m = %$m",
                    normalize = true, 
                    color = 4,
                    xlims = (0, 60), ylims = (0, 0.1),
                    inset = (1, bbox(0.0, 0.0, 0.6, 0.45, :top, :right)),
                    subplot = 2,
                    legend = :topright
                )
                Γ = Γ_calc(J₀, maximum(curr_bonds))
                #println("Γ is $Γ")
                plot!(plt[2], β -> IRFP_density(β, Γ), 0.0:0.1:60.0, label = L"\textrm{IRFP}")
                #println("Added IRFP inset.")
                xlabel!(plt, L"\beta")
                ylabel!(plt, L"P(\beta)")
                savefig(plt, img_path * "IRFP_plot.tikz")
            end
        end

        if verbose
            next!(p)
        end
    end

    # Add final singlet pair to singlets 
    curr_nodes = @view nodes[active]
    singlets[end] = (first(curr_nodes), last(curr_nodes))
end

function sdrg_open!(singlets, system, params)
    #Overwite singlets with output of SDRG process

    # Unpack parameters and data
    @unpack L, plots, boundaries, verbose, img_path, data_path, Δ = params
    @unpack nodes, active, bonds = system

    #Initialise plot
    if plots
        plt = plot(dpi = 300, ylims = (0, 1), xlims = (0, 15), legend= :bottomright)
        hist_points = (round.(Ref(Int64), L/2 .* [0.1, 0.5, 0.7, 0.99]))
    end

    #main loop
    if verbose 
        #println("Bond type is $(eltype(bonds)).")
        p = Progress(L ÷ 2, "IRFP Progress:")
    end

    n_active = sum(active)
    m = 1

    #Check bond provided are odd length 
    if length(bonds) % 2 == 1
        push!(bonds, zero(eltype(bonds)))
    end

    # Make a note of J₀ if we need it for plots
    if plots 
        J₀ = maximum(bonds)
    end
        
    while n_active > 2
        curr_bonds = @view bonds[active]
        curr_L = length(curr_bonds) # should this be n_active? don't think so
        curr_nodes = @view nodes[active]

        # Regardless of boundary conditions, J_m is always Jmax
        Jₘ, Jₘ_position = findmax(curr_bonds)

        # Calculate adjacent spins - useful even in OBC
        Jₗ_position = chain_wrap(Jₘ_position - 1, curr_L)
        Jᵣ_position = chain_wrap(Jₘ_position + 1, curr_L)
        Jₗ = curr_bonds[Jₗ_position]
        Jᵣ = curr_bonds[Jᵣ_position]

        # If J_m is at the start or the end AND we are OBC, need to remove 
        # that pair and move on 

        if Jₘ_position == 1 # Left boundary

            #Deactivate sites at Jₘ, Jᵣ
            curr_active = @view active[active]
            curr_active[Jₘ_position] = false
            curr_active[Jᵣ_position] = false

        #add singlet 
            singlets[m] = (curr_nodes[Jₘ_position], curr_nodes[Jᵣ_position])

        elseif Jₘ_position == curr_L - 1 # right boundary
            #Deactivate sites at Jₘ, Jᵣ
            curr_active = @view active[active]
            curr_active[Jₘ_position] = false
            curr_active[Jᵣ_position] = false
            curr_bonds[Jₗ_position] = zero(eltype(curr_bonds)) #set final bond to zero

            singlets[m] = (curr_nodes[Jₘ_position], curr_nodes[Jᵣ_position])
        else #no boundary 
            J_prime = J′_calc(Jₗ, Jᵣ, Jₘ, Δ)

            # Set Jₗ to J_prime
            curr_bonds[Jₗ_position] = J_prime

            #Deactivate sites at Jₘ, Jᵣ
            curr_active = @view active[active]
            curr_active[Jₘ_position] = false
            curr_active[Jᵣ_position] = false

            #add singlet 
            singlets[m] = (curr_nodes[Jₘ_position], curr_nodes[Jᵣ_position])
        end

        
        m += 1

        n_active = sum(active)

        if plots && in(m, hist_points)
            β = β_calc(maximum(curr_bonds), curr_bonds)
            histogram!(plt, view(β, β .<= 15.0), label = L"m = %$m", 
            normalize = true, bins = 30, xlims = (0, 15))


            if m == last(hist_points)
                #println("Adding IRFP inset...")
                histogram!(plt,
                    β, label = L"m = %$m",
                    normalize = true, 
                    color = 4,
                    xlims = (0, 60), ylims = (0, 0.1),
                    inset = (1, bbox(0.0, 0.0, 0.6, 0.45, :top, :right)),
                    subplot = 2,
                    legend = :topright
                )
                Γ = Γ_calc(J₀, maximum(curr_bonds))
                #println("Γ is $Γ")
                plot!(plt[2], β -> IRFP_density(β, Γ), 0.0:0.1:60.0, label = L"\textrm{IRFP}")
                #println("Added IRFP inset.")
                xlabel!(plt, L"\beta")
                ylabel!(plt, L"P(\beta)")
                savefig(plt, img_path * "IRFP_plot.tikz")
            end
        end

        if verbose
            next!(p)
        end
    end

    # Add final singlet pair to singlets 
    curr_nodes = @view nodes[active]
    singlets[end] = (first(curr_nodes), last(curr_nodes))
end