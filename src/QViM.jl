module QViM

include("Utils.jl")
include("QuantInstr.jl")
include("VirtualMachine.jl")

using .VirtualMachine
export QVM, showstate, step!, execute!

using .QuantInstr
export @quip, QInstr

end


