using Graphs, GraphRecipes

export plot_phase, chord_plot, chord_plot_grid, generate_graph

function generate_graph(singlets)
    L = length(singlets) * 2
    G = SimpleGraph(L, 0)

    for singlet in singlets 
        add_edge!(G, singlet)
    end
    
    G
end

function plot_phase(links; L = length(links) * 2, names = [],
    x = nothing, y = nothing, method = :stress)
    G = generate_graph(links)
    if isnothing(x)
        xs = cos.(range(0, 2π, L + 1))[begin:end - 1]
    else 
        xs = x
    end

    if isnothing(y)
        ys = sin.(range(0, 2π, L + 1))[begin:end - 1]
    else
        ys = y
    end

    graphplot(G, nodesize = 15/L, dpi = 300,
        x = xs, y = ys, names = names, method = method)
end


#=
NOTES
δ = 0.0
=#
function chord_plot(L, δ, h)
    J = randbow_links(δ, h, L)

    ps = SDRGParams(L, false, J)
    out = sdrg(ps)

    G = SimpleGraph(L, 0)
    for link in out.links 
        add_edge!(G, link)
    end 

    graphplot(G, method = :chorddiagram, title = "δ = $δ, h = $h", dpi = 300)
end

function chord_plot_grid(L)
    params = [(0.0, 0.5), (0.0, 10.0), (1.0, 0.5), (1.0, 10.0)]
    plots = Vector{Plots.Plot{Plots.GRBackend}}()
    singlets = fill((0, 0), L ÷ 2)

    for param in params 
        δ, h = param
        distribution = RandbowDistribution()

        ps = SDRGParams(L, false, :open, false, "", "", Δ = 0.0)
        system = generate_system(Js, L)
        out = sdrg(ps)

        G = SimpleGraph(L, 0)
        for link in out.links 
            add_edge!(G, link)
        end 

        push!(plots, graphplot(G, #names = 1:L,
        method = :chorddiagram, title = "δ = $δ, h = $h"))
    end

    plot(plots..., layout = (2, 2), dpi = 300)
end

function saveplot(plt, img_path, name)
    savefig(plt, img_path * name * ".tikz")
    savefig(plt, img_path * "pngs/" * name * ".png")
end

function plotmatrix(M)
    heatmap(Matrix{Float64}(M), yflip = true)
end