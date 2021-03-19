include("QViM.jl")

using .QViM
using LinearAlgebra

quip1 = @quip [
    X >> 1;
    X >> 2 | 1;
]

quip2 = @quip [
    X >> 1;
    X >> 2 | !1;
]

quip3 = @quip [
    X >> 2 | 1;
]

quip4 = @quip [
    1 => C[1]     ;
    X >> 2 | !C[1];
]

quips = [quip1, quip2, quip3, quip4]
for quip in quips
    vm = QVM(quip, 2)
    execute!(vm)
    for (i, j) in enumerate([norm(x)^2 for x in vm.wfn])
        println("$i: $j")
    end
    println()
end
