### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 5bc68900-7ecf-11eb-0d3a-f7f997d8e14a
begin
	using LinearAlgebra
	using Printf
	using PlutoUI
end

# ╔═╡ f3cbdd70-836d-11eb-19d2-25a48489f506
Base.show(io::IO, f::Float64) = @printf io "%.3f" f

# ╔═╡ b90f4c50-7f79-11eb-1552-e997256fd09f
md"
## MapView and it's application to simulating quantum computations
"

# ╔═╡ 79824ea0-7f7b-11eb-2ef3-4d04e143ecfb
md"
#### 1. MapView's declaration and interfaces

MapView is similar to a regular Julia view. However, the data it accesses is determined by an internal map. This allows for a reorganization of the array in question. Used in cases where it is convenient to have two different orderings of the same array at the same time.
"

# ╔═╡ 722db420-7ecf-11eb-063d-b359e7a9020b
struct MapView{T, N} <: AbstractArray{T, N}
	indices::Array{Int, N}
	original::Array{T, N}
end

# ╔═╡ daa83220-7ed2-11eb-0180-bd81ab5dd09d
# AbstractArray interfaces for MapView
begin
	Base.IndexStyle(::Type{MapView}) = IndexLinear()
	Base.getindex(X::MapView, i::Int) = X.original[X.indices[i]]

	function Base.setindex!(X::MapView{T}, v::V, i::Int) where {T, V}
		X.original[X.indices[i]] = v
	end

	Base.firstindex(X::MapView) = 1
	Base.lastindex(X::MapView) = length(X)
end

# ╔═╡ f1909b30-7ed2-11eb-1d69-0371f218c6c1
# Iterator interfaces for MapView
begin
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
end

# ╔═╡ cc9bb1c0-7f7c-11eb-09a9-bb1854c107d8
md"
#### 2. A more efficient way to simulate quantum computation

Below is a use case of MapView, at the same time showcasing a more efficient way to simulate quantum gate operations on a set of qBits.
"

# ╔═╡ 3e690bc0-7eee-11eb-346e-4573f4f5fbd0
⊗(x, y) = kron(x, y)

# ╔═╡ ffd7b6b0-7ef0-11eb-2bdd-67ef92238d66
# Generates a mapping that moves qBit k in N qBits to the top
# If k is negative, treat k as NOT'ed
# O(2^N) complexity
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

# ╔═╡ 69d6e370-7ef5-11eb-3a2b-a7802a5adc8b
# Generates a mapping such that the qbits in ks are at the top of N qBits
# With the top to bottom ordering corresponds to begin to end of ks
# Example:
# 	Let there be 5 qBit lines in order: 
#		ψ₁ ψ₂ ψ₃ ψ₄ ψ₅
#   multibubble([2, 4], 5) would generate a mapping that corresponds to:
#		ψ₂ ψ₄ ψ₁ ψ₃ ψ₅
#
# Complexity: O(n * 2^N + n^2)
# where:
#  		N is the number of qBits
#       n is the number of qBits to be bubbled

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

# ╔═╡ 8d46e6ce-7f99-11eb-1cb4-3358ca6564b5
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

# ╔═╡ 2e692000-7fbd-11eb-0778-d319c4e81810
md"
#### 3. Examples
"

# ╔═╡ 1f485320-7eeb-11eb-2088-1d5554d6d065
# Initialize 3 qBits
begin
	q₁ = [1, 1] / sqrt(2)
	q₂ = [1, 1] / sqrt(2)
	q₃ = [1, -1] / sqrt(2)
end

# ╔═╡ fbfe7840-7fad-11eb-3160-c994a748ae2a
# This routine below simulates the following configuration
#
# |ψ₁⟩ -------x---------
#             |
# |ψ₂⟩ -------|---------
#             |
# |ψ₃⟩ -------x---------
begin
	ψ = q₁ ⊗ q₂ ⊗ q₃
	# This is the swapping operator
	SWAP = [1 0 0 0;
		    0 0 1 0;
		    0 1 0 0;
		    0 0 0 1]
	
	applygate!(ψ, SWAP, [3, 1])
	ψ == q₃ ⊗ q₂ ⊗ q₁
end

# ╔═╡ 45191f7e-7fb3-11eb-1a83-49475fe38f82
# This routine below simulates the following configuration
#
# |ψ₁⟩ -------•---------
#             |
# |ψ₂⟩ -------|---------
#             |
# |ψ₃⟩ -------⊗---------
begin
	ψ2 = q₁ ⊗ q₂ ⊗ q₃
	ψ3 = q₁ ⊗ q₂ ⊗ q₃
	NOT = [0 1;
		   1 0]
	
	applygate!(ψ2, NOT, [3], control=[1])
	applygate!(ψ3, NOT, [3])
	ψ2, ψ3
end

# ╔═╡ 881ecd30-7f85-11eb-3d93-55b1a6dd7039
md"
#### 4. Complexity analysis

##### 4.1. Complexity

Let: 
*  $N$ be the number of qBits
*  $n$ be the size of the quantum gate

First, the generation of a `MapView` of $\psi$ involves the function `multibubble`, which has a complexity of $\mathcal{O}(n2^N + n^2)$

The application of the quantum gate involves a matrix to vector multiplication with a matrix size of $2^{2n}$ and the same number of operations. However, each operation is actually a scalar to vector multiplication with vector size of $2^{N-n}$, which results in $\mathcal{O}(2^{N+n})$ complexity.

Since $2^{N+n} > n2^{N} > n^2$ for large enough $N$ and $n$, the complexity of the whole routine is $\mathcal{O}(2^{N+n})$.

##### 4.2. Comparison to the naive implementation

With the naive implementation, the application of the quantum gate involves matrix to vector multiplications with matrix size $2^{2N}$. This leads to a complexity of $\mathcal{O}(2^{2N})$.

Because $n \leq N$, the above routine is more efficient than the naive implementation. How much of a performance improvement depends on how many qBits involved in an operation compared to the total number of qBits. In the case of $n = N$, that is if the quantum gate involves all available qBits, this method performs no better than the naive implementation. That said, most basic operations only operate on a couple of qBits at a time, and bigger operations can usually be decomposed to a series of smaller ones.

This method is also compatible with a different method that optimizes controled gates, which is not specified in this notebook. This is to be addressed in the future.

"

# ╔═╡ e77cd580-80eb-11eb-2686-016c6fdd8bf1
md"
#### 5. Implementation of the quip syntax

"

# ╔═╡ fa0e3580-8164-11eb-05da-d10765d23a9a
const GATES = Dict(:I => () -> [1 0;
								0 1],

				   :X => () -> [0 1;
								1 0],
				   :RX => γ -> exp(im * γ * GATES[:X]()),


				   :Y => () -> [0 -im;
								im 0],
				   :RY => β -> exp(im * β * GATES[:Y]()),


				   :Z => () -> [1  0;
								0 -1],
				   :RZ => α -> exp(im * α * GATES[:Z]()),


				   :H => () -> [1  1;
								1 -1] / sqrt(2),

				   :S => () -> [1 0;
						        0 im],

				   :T => () -> [1 0;
						        0 exp(im * π / 4)],

				   :PHASE => θ -> [1 0; 0 exp(im * θ)],
				
				   :SWAP => () -> [1 0 0 0;
		                           0 0 1 0;
		 						   0 1 0 0;
		                           0 0 0 1],
)

# ╔═╡ d233fcf0-813e-11eb-0a37-e54e39a36073
# Quantum Instruction each gate represents an operation on the state of the wave function
struct QInstr
	op::Symbol             # Operation code
	func::Function         # Store the subroutine when op == :callfun
	
	param::Vector{Float64} # Parameters for op or func
	
	# Specify qBits that are acted on
	# In measurement mode, specify [target qbit, target cbit]
	target::Vector{Int}    
	
	# specify control bits
	# negative value indicates control-on-zero
	control::Vector{Int}   
	cbit_ctrl::Bool        # True means the control bits are classical 
end

# ╔═╡ 43811d70-813f-11eb-3614-1315573f2d76
# Show interface for gate
function Base.show(io::IO, g::QInstr)
	op = g.op
	param = g.param
	
	# Function call mode
	if op == :funcall
		funcall = Expr(:call, g.func, param...)
		show(io, "$funcall")
	
	# Measurement mode
	elseif op == :measure
		qbit = g.target[1]
		cbit = Expr(:ref, :C, g.target[2])
		show(io, "$qbit => $cbit")
	# Standard gate mode
	else
		op_param = isempty(param) ? op : Expr(:call, op, param...)

		target = length(g.target) == 1 ? g.target[1] : tuple(g.target...)
			
		ctrl_expr = g.cbit_ctrl ? [Expr(:ref, :C, i) for i in g.control] : g.control
		control = length(ctrl_expr) == 1 ? ctrl_expr[1] : tuple(ctrl_expr...)

		if isempty(ctrl_expr)
			show(io, "$op_param >> $target")
		else
			show(io, "$op_param >> $target | $control")
		end
	end
end

# ╔═╡ 61436fd0-7ed7-11eb-23d6-ed68bf61a769
# Show interface for MapView
Base.show(io::IO, x::MapView) = Base.show(io, collect(x))

# ╔═╡ d0eb8520-8143-11eb-1f1b-27c008fd716d
# Functions that helps streamline the parsing procedure
begin
	force_params(s::Symbol) = Expr(:call, s)
	force_params(e::Expr) = e
	
	force_tuple(n::Number) = Expr(:tuple, n)
	force_tuple(s::Symbol) = Expr(:tuple, s)
	function force_tuple(e::Expr)
		if e.head == :tuple
			e
		else
			Expr(:tuple, e)
		end
	end
	
	# Flatten instruction collections
	flatten(i::QInstr) = [i]
	function flatten(a::AbstractArray)
		flattened = [flatten(i) for i in a]
		vcat(flattened...)
	end
	
	dummyfun() = 1 # dummy function
end

# ╔═╡ 7f4adfc0-85aa-11eb-288e-fdf2569ed30b
# Parsing classical bit expressions: C[i]
begin
	function parse_cbit(expr::Expr)
		@assert expr.head == :ref "Not a classical bit: $expr"
		@assert expr.args[1] == :C "Bits can only take the form C[i]"
		expr.args[2]
	end
	
	parse_cbit(other) = throw("Not a classical bit: $other")
end

# ╔═╡ 44603d30-886e-11eb-1313-4d8adc9702da
# Parsing "not" expressions: !Expr
begin
	function parse_not(expr::Expr)
		if expr.head == :call && expr.args[1] == :!
			(-1, expr.args[2])
		else
			(1, expr)
		end
	end
	
	parse_not(other) = (1, other)
end

# ╔═╡ d89eb340-8153-11eb-1341-6518040a962d
begin
	function parse_entry(entry::Expr)
		##############################################
		# List expresion:
		# 	[sub expression 1, sub expression 2, ...]
		if entry.head == :vcat || entry.head == :vect
			parsed = [parse_entry(e) for e in entry.args]
			return Expr(:vcat, parsed...)
		
		##############################################
		# List comprehension expression:
		# 	[sub_expr for_expr if_expr]
		elseif entry.head == :comprehension
			generator = entry.args[1]
			generator.args[1] = parse_entry(generator.args[1])
			generator.args[2] = esc(generator.args[2])
			return Expr(:comprehension, generator)
		
		###############################
		# QInstr expression 
		# 	+ Measure mode
		#	+ Standard gate mode
		#	+ Function call mode
		elseif entry.head == :call
			control = []
			classical = []
			cbit_ctrl = false
			funcall = dummyfun # dummy function
			
			####################################
			# Measure mode: qBit => cBit
			if entry.args[1] == :(=>)
				op = Meta.quot(:measure)
				qbit = entry.args[2]
				cbit = parse_cbit(entry.args[3])
				target = Expr(:vect, qbit, cbit)
				
				return :(QInstr($op, $funcall, [], $target, $control, $cbit_ctrl))
			end
			
			###################################################
			# Standard gate mode: OP(param) >> target | control
			if entry.args[1] == :| # Parsing control
				control_expr = force_tuple(entry.args[3]).args
				
				# Parsing NOT'ed expressions
				control_expr = parse_not.(control_expr)
				
				# Check for classical control mode
				if typeof(control_expr[1][2]) == Expr
					cbit_ctrl = control_expr[1][2].head == :ref
				end
				
				# Process control bit expressions depends on cbit_mode
				if cbit_ctrl
					# Parse cBit values
					control_expr = [(mult, parse_cbit(expr)) 
						for (mult, expr) in control_expr]
					
				# qbit mode
				else
					# sanity check
					for expr in control_expr
						if typeof(expr) == Expr
							@assert expr.head != :ref "not a qBit $expr"
						end
					end
				end
				# Applying NOT'ed expressions
				control_expr = [Expr(:call, :*, mult, expr)
					for (mult, expr) in control_expr]
				
				control = esc(Expr(:vect, control_expr...))
				entry = entry.args[2]
			end
			
			@assert entry.head == :call "Invalid quip entry $entry"
			if entry.args[1] == :>> # Parsing Operation
				# Parse target
				target = esc(Expr(:vect, force_tuple(entry.args[3]).args...))

				# Parse op and param
				op_and_param = force_params(entry.args[2]).args
				op = Meta.quot(op_and_param[1])
				param = esc(Expr(:vect, op_and_param[2:end]...))
				
				return :(QInstr($op, $funcall, $param, $target, $control, $cbit_ctrl))
			end
			
			##################################
			# Function call mode: func(param)
			op_and_param = entry.args
			op = Meta.quot(:funcall)
			param = esc(Expr(:vect, op_and_param[2:end]...))
			
			func = esc(op_and_param[1])
			
			return :(QInstr($op, $func, $param, [], $control, $cbit_ctrl))

		# Branch expression a ? b : c
		elseif entry.head == :if
			a = entry.args[1]
			b = entry.args[2]
			c = entry.args[3]
			return Expr(:if, esc(a), parse_entry(b), parse_entry(c))
			
		# Dot expression a.b
		elseif entry.head == :(.)
			return entry
		end
		
		throw("can't parse $entry")
	end
	
	function parse_entry(entry::Symbol)
		esc(entry)
	end
end

# ╔═╡ 317c8920-813f-11eb-2048-f1ff2bafeab5
macro quip(expr)
	if expr.head == :(=)
		lhs = expr.args[1]
		rhs = expr.args[2]
	else
		rhs = expr
	end
	
	if rhs.head == :block
		rhs = rhs.args[2]
	end
	
	rhs = parse_entry(rhs)

	rhs = Expr(:call, :flatten, rhs)
	if expr.head == :(=)
		:($(esc(lhs)) = $rhs)
	else
		:($rhs)
	end
end

# ╔═╡ 5a916b00-80a9-11eb-2cc6-93f42781e1b6
####################################
# Quip declaration example
begin
	quip1(rx_param) = @quip [
		RX(rx_param) >> 3 | 1  ;
		H            >> 1 | 2  ;
		1 => C[2]]
	
	quip2 = @quip [
	 	H >> 2;
	 	quip1(1)]

	quip3(a) = @quip [iseven(i) ? H >> (i % 3 + 1) : quip2 for i in 1:a]
	
	(quip1(π / 4), quip2, quip3(10))
end

# ╔═╡ c73aafc0-8287-11eb-1a3e-2de3296a35dc
# Type aliases
const Quip = Vector{QInstr}

# ╔═╡ 114b56c0-8240-11eb-15bc-dd322c881361
# The Quantum Virtual Machine, 
# 	a Turing machine that holds the instruction set
#	as well as the current state of the wave function
struct QVM{R <: Real}
	# Stacks, each element is a stack frame
	stack::Vector{Quip}
	pc_stack::Vector{Int}
	
	# Classical memory
	reg::Vector{Bool}
	
	# Wavefunction
	wfn::Vector{Complex{R}}
	 
	# A dictionary of primitive operations
	dict::Dict{Symbol, Function}
	
	function QVM{R}(quip::Vector{QInstr}, nbit::Int; 
			dict=GATES, nreg=32 ) where R <: Real
		
		wfn = zeros(Complex{R}, 2^nbit)
		wfn[1] = 1
		
		registers = zeros(Int, nreg)
		new([quip], [1], registers, wfn, dict)
	end
end

# ╔═╡ 6929cfb0-82c8-11eb-3857-8df1ea9a3a88
# pretty print the current stack
# for internal use by showstate
function showstack(stack::Vector{QInstr}, pc::Int, n::Int; io::IO=stdout)
	len = length(stack)
	
	# Checking boundaries
	if pc - n < 1
		s = 1
		e = 1 + 2n > len ? len : 1 + 2n
	elseif pc + n > len
		e = len
		s = 1 - 2n < 1 ? 1 : len - 2n
	else
		s = pc - n
		e = pc + n
	end
	
	# Print stack
	if s != 1
		println(io, "...")
	end
	
	for i in s:e
		print(io, "$i: $(stack[i])")
		if i == pc
			println(io, "   <--- pc")
		else
			println(io)
		end
	end
	
	if e != len
		println(io, "...")
	end
	
	if pc > len
		println(io, "$(len + 1):  <- pc")
	end
end

# ╔═╡ 4faadd40-82c8-11eb-2d7e-f7e533a7e495

function showstate(vm::QVM, n=2; io::IO=stdout)
	println(io, "##########################")
	for i in eachindex(vm.stack)
		stack = vm.stack[i]
		pc = vm.pc_stack[i]
		
		println(io, "$i ---------------------------")
		showstack(stack, pc, n, io=io)
		println(io)
	end
end

# ╔═╡ f593d63e-823b-11eb-2f6b-a91be5d313bf
function step!(vm::QVM)
	dict = vm.dict
	quip = vm.stack[end]
	pc = vm.pc_stack[end]
	
	################################################
	# PC reached the end of the current stack frame
	#
	if pc > length(quip)
		# VM reached the end of program do nothing
		if length(vm.stack) == 1
			return false
		end
		
		# pop stack frame
		pop!(vm.stack)
		pop!(vm.pc_stack)
		vm.pc_stack[end] += 1 # increment pc
		return true
	end

	instr = quip[pc]
	op = instr.op
	func = instr.func
	param = instr.param
	target = instr.target
	control = instr.control
	cbit_ctrl = instr.cbit_ctrl
	
	
	#####################################
	# Function call push new stack frame
	#
	if op == :funcall
		new_quip = func(param...)
		push!(vm.stack, new_quip)
		push!(vm.pc_stack, 1)
	
	##########################################################
	# Measurement, collapse wave function and save to register
	#
	elseif op == :measure
		qbit = target[1]
		cbit = target[2]
		# Bases
		bas0 = [1; 0]
		bas1 = [0; 1]
		
		wfn0 = deepcopy(vm.wfn)
		applygate!(wfn0, bas0 * bas0', [qbit])
		P0 = vm.wfn' * wfn0
		
		# collapse to 0 state
		if rand() < norm(P0)
			vm.reg[cbit] = false
			vm.wfn .= wfn0 / sqrt(P0)

		# collapse to 1 state
		else
			wfn1 = vm.wfn - wfn0
			P1 = 1 - P0
			
			vm.reg[cbit] = true
			vm.wfn .= wfn1 / sqrt(P1)
		end
		vm.pc_stack[end] += 1 # increment pc
	
	##################################
	# Gate operation, call applygate!
	#
	else
		operator = get(dict, op) do
			throw("$op not in dict, this should not happen")
		end
		
		# cBit control
		if cbit_ctrl
			conds = [vm.reg[abs(i)] ⊻ (i < 0) for i in control]
			if reduce((a, b) -> a && b, conds)
				applygate!(vm.wfn, operator(param...), target)
			end
		
		# qBit control
		else
			applygate!(vm.wfn, operator(param...), target, control=control)
		end
		vm.pc_stack[end] += 1 # increment pc
	end
		
	
	true
end
	

# ╔═╡ ea7a2ed0-8600-11eb-2497-bb23e393d273
begin
	qft_rotations(n, k, o) = @quip [j==k ? 
										H >> k : 
										PHASE(2π / 2^(j-k+1)) >> k | j 
									for j in k:(n+o)]
	
	qft_swaps(n, o) = @quip [SWAP >> (j+o, n-j+1+o) for j in 1:div(n, 2)]
	
	qft(n, o=0) = @quip [qft_rotations(n, k, o) for k in (1+o):(n+o)]
	
end

# ╔═╡ a59e1620-828e-11eb-2a19-15b99805ac4d
with_terminal() do
	N = 2
	quip = @quip [
		X >> 1            ;
		1 => C[1]         ;
		X >> 2 | !C[1]    ;
	]
	
	for i in 1:1
		vm = QVM{Float64}(quip, N)
		while step!(vm)
			# showstate(vm)
		end
		for (i, j) in enumerate([norm(x)^2 for x in vm.wfn])
			println("$i: $j")
		end
	end
end

# ╔═╡ ae503700-88db-11eb-2c85-5574fe9aeca5
 @quip [i == 3 ? H >> i : X >> i | (i + 1) for i in 1:3]

# ╔═╡ Cell order:
# ╠═5bc68900-7ecf-11eb-0d3a-f7f997d8e14a
# ╠═f3cbdd70-836d-11eb-19d2-25a48489f506
# ╟─b90f4c50-7f79-11eb-1552-e997256fd09f
# ╟─79824ea0-7f7b-11eb-2ef3-4d04e143ecfb
# ╠═722db420-7ecf-11eb-063d-b359e7a9020b
# ╠═daa83220-7ed2-11eb-0180-bd81ab5dd09d
# ╠═f1909b30-7ed2-11eb-1d69-0371f218c6c1
# ╠═61436fd0-7ed7-11eb-23d6-ed68bf61a769
# ╟─cc9bb1c0-7f7c-11eb-09a9-bb1854c107d8
# ╠═3e690bc0-7eee-11eb-346e-4573f4f5fbd0
# ╠═ffd7b6b0-7ef0-11eb-2bdd-67ef92238d66
# ╠═69d6e370-7ef5-11eb-3a2b-a7802a5adc8b
# ╠═8d46e6ce-7f99-11eb-1cb4-3358ca6564b5
# ╟─2e692000-7fbd-11eb-0778-d319c4e81810
# ╠═1f485320-7eeb-11eb-2088-1d5554d6d065
# ╠═fbfe7840-7fad-11eb-3160-c994a748ae2a
# ╠═45191f7e-7fb3-11eb-1a83-49475fe38f82
# ╟─881ecd30-7f85-11eb-3d93-55b1a6dd7039
# ╟─e77cd580-80eb-11eb-2686-016c6fdd8bf1
# ╠═fa0e3580-8164-11eb-05da-d10765d23a9a
# ╠═d233fcf0-813e-11eb-0a37-e54e39a36073
# ╠═43811d70-813f-11eb-3614-1315573f2d76
# ╠═d0eb8520-8143-11eb-1f1b-27c008fd716d
# ╠═7f4adfc0-85aa-11eb-288e-fdf2569ed30b
# ╠═44603d30-886e-11eb-1313-4d8adc9702da
# ╠═d89eb340-8153-11eb-1341-6518040a962d
# ╠═317c8920-813f-11eb-2048-f1ff2bafeab5
# ╠═5a916b00-80a9-11eb-2cc6-93f42781e1b6
# ╠═c73aafc0-8287-11eb-1a3e-2de3296a35dc
# ╠═114b56c0-8240-11eb-15bc-dd322c881361
# ╠═6929cfb0-82c8-11eb-3857-8df1ea9a3a88
# ╠═4faadd40-82c8-11eb-2d7e-f7e533a7e495
# ╠═f593d63e-823b-11eb-2f6b-a91be5d313bf
# ╠═ea7a2ed0-8600-11eb-2497-bb23e393d273
# ╠═a59e1620-828e-11eb-2a19-15b99805ac4d
# ╠═ae503700-88db-11eb-2c85-5574fe9aeca5
