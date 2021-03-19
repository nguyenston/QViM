include("QViM.jl")

using .QViM

quip = @quip [
    H >> 3;
    H >> 2;
    H >> 1;
]

println(quip)
