module QuantInstr

using ..Utils

export @quip, QInstr

"""
Quantum Instruction each gate represents an operation on the state of the wave function
"""
struct QInstr
	op::Symbol             # Operation code
	func::Function         # Store the subroutine when op == :callfun
	
	param::Vector{Float64} # Parameters for op or func
	
	# Specify qBits that are acted on
	# In measurement mode, specify [target qbit, target cbit]
	target::Vector{Int}    
	
	control::Vector{Int}   # specify control bits
	cbit_ctrl::Bool        # True means the control bits are classical 
end

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

		if isempty(control)
			show(io, "$op_param >> $target")
		else
			show(io, "$op_param >> $target | $control")
		end
	end
end

#######################################################
# Functions that helps streamline the parsing procedure
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


#########################################
# Parsing classical bit expressions: C[i]
function parse_cbit(expr::Expr)
    @assert expr.head == :ref "Not a classical bit: $expr"
    @assert expr.args[1] == :C "Bits can only take the form C[i]"
    return expr.args[2]
end

parse_cbit(other) = throw("Not a classical bit: $other")

"""
Parses a QInsr expression
"""
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
            
            # Check for classical control mode
            if typeof(control_expr[1]) == Expr
                cbit_ctrl = control_expr[1].head == :ref
            end
            
            # Process control bit expressions depends on cbit_mode
            if cbit_ctrl
                control_expr = [parse_cbit(expr) for expr in control_expr]
            else
                for expr in control_expr
                    if typeof(expr) == Expr
                        @assert expr.head != :ref "not a qBit $expr"
                    end
                end
            end
            control = esc(Expr(:vect, control_expr...))
            entry = entry.args[2]
        end
        
        @assert entry.head == :call "Invalid quip entry $entry"
        if entry.args[1] == :>> # Parsing Operation
            target = esc(Expr(:vect, force_tuple(entry.args[3]).args...))

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
    end
    
    throw("can't parse $entry")
end

function parse_entry(entry::Symbol)
    esc(entry)
end

"""
Usage:
@quip <quip expression>
"""
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

end
