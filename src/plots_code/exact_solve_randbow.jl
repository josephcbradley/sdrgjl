using ArnoldiMethod, SparseArrays, BandedMatrices, Polyester, DataFrames, ProgressMeter, OnlineStats

struct ExactBuffer{T}
    J::Vector{T}
    T::SparseMatrixCSC{T, Int64}
    C::Matrix{T}

    ExactBuffer(L::Int, data_type::DataType) = new{data_type}(
        ones(data_type, L),
        SparseMatrixCSC(BandedMatrix(-1 => ones(data_type, L - 1), 1 => ones(data_type, L - 1))),
        Matrix{data_type}(undef, L, L)
    )
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

L = 100
L_over_2 = L ÷ 2
setprecision(BigFloat, 128)
data_type = BigFloat
l_range = 2:1:32
ps = [(1.5, 0.5), (3, 1), (4.5, 1.5)]
plots = Plots.Plot{Plots.PGFPlotsXBackend}[]
plt = plot(dpi = 600, xlabel = L"l", ylabel = L"S", legend = :bottomright,
            xscale = :log2, yscale = :log2)

exact_trials = 1000

for parameters in reverse(ps)
    h, d = parameters

    distribution = RandbowDistribution{data_type}(d, h)
    ls = collect(l_range)
    entropy = n_means(length(l_range))
    #data = DataFrame(l = collect(l_range), entropy = zeros(data_type, length(l_range)))
    p = Progress(exact_trials, 1)
    @batch threadlocal = ExactBuffer(L, data_type)::ExactBuffer{data_type} for trials in 1:exact_trials
            
        J = threadlocal.J
        N = length(J)
        T = threadlocal.T
        C = threadlocal.C

        
        J[1:N-1] .= rand(distribution, L - 1)
        J[N] = zero(data_type)

        #overwrite lower diagonal of T 
        for i in 2:L
            j = i - 1
            T[i, j] = J[j] / 2
        end
        #overwrite upper diagonal of T 
        for i in 1:L-1
            j = i + 1
            T[i, j] = J[i] / 2
        end

        
        #Eigensolve T 
        decomp, history = ArnoldiMethod.partialschur(T, nev = L ÷ 2, which = SR())

        
        # Label eigenvectors of groundstate
        gs_eigvecs = decomp.Q

        #Check eigenvectors
        for i in 1:L_over_2
            v = gs_eigvecs[:, i]
            λ = decomp.R[i, i]
            check = isapprox(norm((T - λ*I(L))*v), 0, atol = 1e-20)
            if !check 
                @warn "Eigenvalues are not accurate!"
            end
        end


        # Calculate values of C
        for  j in 1:L      
            for i in 1:L
                C[i, j] = cij(i, j, gs_eigvecs)
            end
        end


        C .= Symmetric(C)
        
        
        #Now calculate the entropy for all 64 values of l 
        local_entropy = zeros(length(l_range))
        for i in eachindex(ls)
            l = ls[i]
            sub_C = @view C[L_over_2 + 1:L_over_2 + l, L_over_2 + 1:L_over_2 + l]
            sub_C_decomp, sub_C_history = partialschur(sub_C, nev = l)
            sub_C_vals = diag(sub_C_decomp.R)
            #Check eigenvectors
            for i in 1:l
                v = sub_C_decomp.Q[:, i]
                λ = sub_C_vals[i]
                check = isapprox((sub_C - λ*I(l))*v, zeros(data_type, l), atol = 1e-20)
                if !check 
                    @warn "Eigenvalues are not accurate!"
                end
            end
            filter!(x -> x > zero(typeof(x)), sub_C_vals)
            #filter!(x -> !isapprox(x, 0, atol = 1e-2), sub_C_vals)
            local_entropy[i] = eig_entropy(sub_C_vals)

        end
        fit!(entropy, local_entropy)
    
        next!(p) 
    end

    data = DataFrame(l = ls, entropy = Vector{Float64}(undef, length(ls)), copycols = false)
    for i in eachindex(data.entropy)
        data.entropy[i] = value(entropy[i])
    end

    @df data scatter!(plt, :l, :entropy, 
        label = L"h = %$h , \delta = %$d")

    #push!(plots, subplot)
end

saveplot(plt, img_path, "exact_randbow_solve")