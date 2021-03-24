module QViM

include("Utils/Utils.jl")
include("QuantInstr.jl")
include("VirtualMachine.jl")

using .VirtualMachine
export QVM, showstate, step!, execute!

using .QuantInstr
export @quip, QInstr

end


