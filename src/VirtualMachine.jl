module VirtualMachine

export QVM, showstate, step!, execute!

using ..Utils
using ..QuantInstr

⊗(x, y) = kron(x, y)

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



# Type aliases
const QIS = Vector{QInstr}

"""
The Quantum Virtual Machine, 
	a Turing machine that holds the instruction set
	as well as the current state of the wave function
"""
struct QVM	
    # Stacks, each element is a stack frame
	stack::Vector{QIS}
	pc_stack::Vector{Int}
	
	# Classical memory
	reg::Vector{Bool}
	
	# Wavefunction
	wfn::Vector{Complex{Float64}}
	 
	# A dictionary of primitive operations
	dict::Dict{Symbol, Function}
	
	function QVM(qis::Vector{QInstr}, nbit::Int; 
                    dict=GATES, nreg=32 )
		
		wfn = zeros(Complex{Float64}, 2^nbit)
		wfn[1] = 1
		
		registers = zeros(Int, nreg)
		new([qis], [1], registers, wfn, dict)
	end
end

"""
pretty print the current stack
for internal use by showstate
"""
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

function step!(vm::QVM)
	dict = vm.dict
	qis = vm.stack[end]
	pc = vm.pc_stack[end]
	
	################################################
	# PC reached the end of the current stack frame
	#
	if pc > length(qis)
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

	instr = qis[pc]
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
		new_qis = func(param...)
		push!(vm.stack, new_qis)
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
			conds = [vm.reg[i] for i in control]
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

"""
Run step!(vm) until pc reaches the end of the program strip
"""
function execute!(vm::QVM; verbose::Bool=false, io::IO=stdout)
    if verbose 
        showstate(vm, io=io) 
    end
    while step!(vm)
        if verbose 
            showstate(vm, io=io) 
        end
    end
end
	
end
