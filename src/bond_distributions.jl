# Functions for generating particular bond distributions 

export PaolaDistribution, inverse_K_distribution, randbow_links, RandbowDistribution, FixedDistribution

using Distributions

import Distributions.Sampleable, Distributions.Univariate, Distributions.Continuous

struct PaolaDistribution <: Sampleable{Univariate, Continuous}
    δ::Float64
end

PaolaDistribution() = PaolaDistribution(0.1)

import Base.rand

import Random.Random.AbstractRNG

function inverse_K_distribution(r, δ)
    #For r ∈ {0, 1}, return the value of K associated with that 
    # cumulative density r 
    r ^ δ
end

function rand(rng::AbstractRNG, s::PaolaDistribution)
    p = rand(Uniform(0.0, 1.0))
    inverse_K_distribution(p, s.δ)
end

using OffsetArrays

function randbow_J(K, m, h)
    if m != 0
        K * exp(-h * abs(m))
    else
        K * exp(-h / 2)
    end
end

function randbow_links(δ, h, L::T) where T <: Integer
    # Clean case is δ = 0.0
    
    base_out = Vector{typeof(δ)}(undef, L)
    #Create offset wrapper for easy indexing. Can return
    #base out as modifications to offset flow up to parent
    #Two cases - L is even or odd 
    # Easy if L is odd:
    if L % 2 == 1
        J = OffsetArray(base_out, -L ÷ 2:L ÷ 2)
    else
        #if L is even, let the RHS of the chain be longer
        J = OffsetArray(base_out, -L ÷ 2 + 1:L ÷ 2)
    end


    #Produce randbow links for disorder parameter δ
    # and decay parameter h for chain of length L.

    #From Paola 2018:
    # Jₘ = Kₘ x {exp(-h/2) if m == 0, else exp(-h*abs(m))}

    # Prepare vector of Kₘ
    P = PaolaDistribution(δ)
    base_K = rand(P, L)
    # want offset K too 

    if L % 2 == 1
        K = OffsetArray(base_K, -L ÷ 2:L ÷ 2)
    else
        #if L is even, let the RHS of the chain be longer
        K = OffsetArray(base_K, -L ÷ 2 + 1:L ÷ 2)
    end

    for m in eachindex(J)
        if m == 0
            J[m] = K[m] * exp(-h / 2.0)
        else 
            J[m] = K[m] * exp(-h * abs(m))
        end
    end

    base_out
end

struct RandbowDistribution{T}
    δ::T
    h::T
end

function rand(r::RandbowDistribution, N::Integer)
    δ, h = r.δ, r.h
    randbow_links(δ, h, N)
end

struct FixedDistribution 
    data::Vector{Float64}
    L::Int64
end

function FixedDistribution(v)
    L = length(v)
    FixedDistribution(v, L)
end

function rand(r::FixedDistribution, N::Integer)
    if N != r.L
        error("Distribution length ($(r.L)) is different to N ($(N))!")
    end

    return r.data
end

function power_law(i, α)
    if i != 0
        inv(abs(i)^α)
    else
        inv(1^α) * 2
    end
end

function powerbow_links(δ, α, L::T) where T <: Integer
    # Clean case is δ = 0.0
    
    base_out = Vector{typeof(δ)}(undef, L)
    #Create offset wrapper for easy indexing. Can return
    #base out as modifications to offset flow up to parent
    #Two cases - L is even or odd 
    # Easy if L is odd:
    if L % 2 == 1
        J = OffsetArray(base_out, -L ÷ 2:L ÷ 2)
    else
        #if L is even, let the RHS of the chain be longer
        J = OffsetArray(base_out, -L ÷ 2 + 1:L ÷ 2)
    end


    #Produce randbow links for disorder parameter δ
    # and decay parameter h for chain of length L.

    #From Paola 2018:
    # Jₘ = Kₘ x {α / 2 if m == 0, else inv(abs(i)^α) / 2}

    # Prepare vector of Kₘ
    P = PaolaDistribution(δ)
    base_K = rand(P, L)
    # want offset K too 

    if L % 2 == 1
        K = OffsetArray(base_K, -L ÷ 2:L ÷ 2)
    else
        #if L is even, let the RHS of the chain be longer
        K = OffsetArray(base_K, -L ÷ 2 + 1:L ÷ 2)
    end

    for m in eachindex(J)
        J[m] = K[m] * power_law(m, α)
    end

    base_out
end

struct PowerbowDistribution{T<:Real}
    δ::T
    α::T
end

function rand(r::PowerbowDistribution, N::Integer)
    δ, α = r.δ, r.α
    powerbow_links(δ, α, N)
end