struct WrappedVec{T} <: AbstractVector{T}
    data::Vector{T}
end

import Base.resize!, Base.deleteat!    

import Base.getindex, Base.setindex!, Base.setindex, Base.size

function Base.getindex(V::WrappedVec, i::Core.Integer)::eltype(V.data)
    N = length(V)
    if i < 1 
        @inbounds Base.getindex(V.data, N + i)
    elseif i > N
        @inbounds Base.getindex(V.data, i - N)
    else 
        @inbounds Base.getindex(V.data, i)
    end
end

function Base.setindex(V::WrappedVec, X, i::Core.Integer)
    N = length(V)
    if i < 1 
        @inbounds Base.setindex(V.data, X, N + i)
    elseif i > N
        @inbounds Base.setindex(V.data, X, i - N)
    else 
        @inbounds Base.setindex(V.data, X, i)
    end
end

function Base.setindex!(V::WrappedVec, X, i::Core.Integer)
    N = length(V)
    if i < 1 
        @inbounds Base.setindex!(V.data, X, N + i)
    elseif i > N
        @inbounds Base.setindex!(V.data, X, i - N)
    else 
        @inbounds Base.setindex!(V.data, X, i)
    end
end

function Base.size(V::WrappedVec)
    Base.size(V.data)
end    

function Base.deleteat!(w::WrappedVec, i::Core.Integer)
    N = length(w)
    if i < 1 
        @inbounds Base.deleteat!(w.data, N + i)
    elseif i > N
        @inbounds Base.deleteat!(w.data, i - N)
    else 
        @inbounds Base.deleteat!(w.data, i)
    end
end


Base.minimum(w::WrappedVec) = Base.minimum(w.data)
Base.maximum(w::WrappedVec) = Base.maximum(w.data)
Base.argmax(w::WrappedVec) = Base.argmax(w.data)
Base.argmin(w::WrappedVec) = Base.argmin(w.data)

Base.resize!(w::WrappedVec, i::Core.Integer) = Base.resize!(w.data, i)

import Base.filter!

function Base.filter!(f, w::WrappedVec)
    Base.filter!(f, w.data)
end

import Base.eltype

Base.eltype(w::WrappedVec) = Base.eltype(w.data)


function wrap_zero(x, L)
    if x >= zero(eltype(x))
        x
    else
        L + x
    end
end