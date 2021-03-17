# QViM: Quantum Virtual Machine
As the name suggests this package allows for the creation of simulated quantum machines.

## Basic Usage

A quantum circuit can be thought of as a series of gate applications.
Similarly, the virtual machine class `QVM` is essentially a turing machine with an internal state (wave function) that takes in
a list of `QInstr` (named Quantum Program or QuIP) and execute on it procedually.

```julia
using QViM

# Creating a instruction strip
program = @quip [
  H >> 1     ;
  H >> 2 | 1 ;
  X >> 3 | 2 ;
]

# Creating a virtual machine running program with 3 qBits
qvm = QVM(program, 3)

# Run the virtual machine
execute!(qvm)
```
