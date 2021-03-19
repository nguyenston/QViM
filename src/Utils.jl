module Utils

export MapView, applygate!

using LinearAlgebra

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





"""
Generates a mapping that moves qBit k in N qBits to the top
To be used internally by multibubble
O(2^N) complexity
"""
function bubble(k::Int, N::Int)
	index_list = collect(1:2^N)
	part_size = 2^(N - abs(k))
	partitions = collect(Iterators.partition(index_list, part_size))
	
	odds  = 1:2:length(partitions)
	evens = 2:2:length(partitions)
	
	if k > 0
		part_map = collect(Iterators.flatten((odds, evens)))
	# if k < 0 generate a map correspond to the qBit in question being NOT'ed
	else
		part_map = collect(Iterators.flatten((evens, odds)))
	end
	
	mapped_parts = MapView(part_map, partitions)
	collect(Iterators.flatten(mapped_parts))
end

"""
Generates a mapping such that the qbits in ks are at the top of N qBits
With the top to bottom ordering corresponds to begin to end of ks
Negative k indicate control-on-zero
Example:
	Let there be 5 qBit lines in order: 
		ψ₁ ψ₂ ψ₃ ψ₄ ψ₅
  multibubble([2, 4], 5) would generate a mapping that corresponds to:
		ψ₂ ψ₄ ψ₁ ψ₃ ψ₅

Complexity: O(n * 2^N + n^2)
where:
 		N is the number of qBits
      	n is the number of qBits to be bubbled
"""
function multibubble(ks::Vector{Int}, N::Int)
	final_map = collect(1:2^N)
	
	push_down_list = []
	for k in reverse(ks)
		# separate sign from value
		(sign, k) = (div(k, abs(k)), abs(k))
		
		push_down = 0
		for i in push_down_list
			if k < i
				push_down += 1
			elseif k == i
				throw("bubbling list has dublicate entry")
			end
		end
		push!(push_down_list, k)
			
		map = bubble(sign * (k + push_down), N)
		final_map = collect(MapView(map, final_map))
	end
	
	final_map
end

"""
TODO: document
"""
function applygate!(
		ψ::V, 
		Gate::M, 
		target::Vector{Int};
		control::Vector{Int} = Vector{Int}()
		) where {V <: AbstractVector{<: Number}, M <: AbstractMatrix{<: Number}}
	
	N = Int(log2(length(ψ)))
	n = length(target)
	
	@assert size(Gate)[1] == size(Gate)[2] "Invalid gate"
	@assert size(Gate)[1] == 2^n "Mismatch gate and # of input"
	
	top_qbit_list = [control; target]
	partition_size = div(length(ψ), 2^length(top_qbit_list))

    @assert max(top_qbit_list...) <= N "qBit number larger than number of qBits"
	
	# Generate a mapping that bubbles all relevent qBits to the top
	# Since in controlled gates' operations,
	# 	only the bottom right of the matrix (the Gate application) is non-trivial,
	#   only takes the end side of ψ that corresponds to the Gate 
	bubble_map = multibubble(top_qbit_list, N)
	gate_index = length(ψ) - size(Gate)[1] * partition_size + 1
	map_ψ = MapView(bubble_map[gate_index:end], ψ)
	
	# Apply Gate on ψ
	map_ψ_part = collect(Iterators.partition(map_ψ, partition_size))
	map_ψ .= collect(Iterators.flatten(Gate * map_ψ_part)) # ψ mutated here
end

end
