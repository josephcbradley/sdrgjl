using Parameters

struct System{R, T}
    nodes::R
    active::Vector{Bool}
    bonds::Vector{T}
end

function generate_system(Js, L)
    # Generate a system given pre-made Js
    System(
        1:L, 
        fill(true, L),
        Js
    )
end

import Base.minimum, Base.maximum, Base.argmax, Base.argmin

Base.minimum(s::System) = Base.minimum(s.bonds)
Base.maximum(s::System) = Base.maximum(s.bonds)
Base.argmax(s::System) = Base.argmax(s.bonds)
Base.argmin(s::System) = Base.argmin(s.bonds)


@with_kw struct SDRGParams{T<:Real}
    L::Int64 = 1000
    plots::Bool = false
    boundaries::Symbol = :periodic
    verbose::Bool = false
    img_path::String = ""
    data_path::String = ""
    Δ::T = 0.0 #Default to XX case 

    #= function SDRGParams(L, plots, bounds, verbose, img_path, data_path)
        if L % 2 != 0
            error("L is odd! L must be even.")
        end

        if bounds ∉ [:periodic, :open]
            error("Your boundary conditions are not recognised!")
        end

        new(L, plots, bounds, verbose, img_path, data_path)
    end =#
end

#= struct DisorderAverageBuffer{T, R}
    Js::Vector{T}
    singlets::Vector{Tuple{Int64, Int64}}
    system::System{R, T}
    buffer::DataFrame
end

function new_disorder_buffer(T, L, variable_name, 
    measurement_name,
    buffer_length)
    DisorderAverageBuffer(
        Vector{T}(undef, L),
        fill((0, 0), L ÷ 2), 
        generate_system(Vector{T}(undef, L), L), 
        DataFrame(Symbol(variable_name) => Vector{Int64}(undef, buffer_length),
                Symbol(measurement_name) => Vector{Float64}(undef, buffer_length)
        )
    )    
end =#

using SparseArrays

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