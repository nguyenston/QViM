module MapViews

using LinearAlgebra

export MapView

"""
MapView is similar to a regular Julia view. However, the data it accesses 
is determined by an internal map. This allows for a reorganization 
of the array in question. Used in cases where it is convenient to have 
two different orderings of the same array at the same time.
"""
struct MapView{T, N} <: AbstractArray{T, N}
    indices::Array{Int, N}
    original::Array{T, N}
end

######################################
# AbstractArray interfaces for MapView
Base.IndexStyle(::Type{MapView}) = IndexLinear()
Base.getindex(X::MapView, i::Int) = X.original[X.indices[i]]

function Base.setindex!(X::MapView{T}, v::V, i::Int) where {T, V}
    X.original[X.indices[i]] = v
end

Base.firstindex(X::MapView) = 1
Base.lastindex(X::MapView) = length(X)

######################################
# Iterator interfaces for MapView
function Base.iterate(iter::MapView)
    if length(iter) == 0
        nothing
    else
        iter[1], 1
    end
end

function Base.iterate(iter::MapView, state)
    if state >= length(iter)
        nothing
    else
        iter[state + 1], state + 1
    end
end

Base.eltype(::Type{MapView{T}}) where T = T
Base.size(iter::MapView) = size(iter.indices)

##############################
# Show interface for MapView
Base.show(io::IO, x::MapView) = Base.show(io, collect(x))

end