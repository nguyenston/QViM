# QViM: Quantum Virtual Machine
As the name suggests this package allows for the creation of simulated quantum machines.

## 1. Basic Usage

A quantum circuit can be thought of as a series of gate applications.
Similarly, the virtual machine class `QVM` is essentially a turing machine with an internal state (wave function) that takes in
a list of `QInstr` (named Quantum Program or QuIP) and execute on it procedually.

```julia
using QViM

# Creating a instruction strip using the @quip macro
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
## 2. QInstr syntax

### 2.1. Basic expression
The expression below (reads _**apply `U` on `a` given `b`**_) represents the following circuit:
```julia
@quip [U >> a | b]
```

<img src="https://github.com/npnguyen99/QViM/blob/main/assets/Basic_QInstr.png?raw=true" height=150>

In the case of multiple target or control qBits, use comma separated tuples.
```julia
@quip [F >> (a1, a2) | (b1, b2)]
```

<img src="https://github.com/npnguyen99/QViM/blob/main/assets/Multibit_QInstr.png?raw=true" height=300>

### 2.2. Function call
In cases where a pattern are repeatedly used, a function call that returns a QuIP can be used as an expression.
```Julia
subroutine(a, b, c) = @quip [
  H >> a     ;
  X >> b | a ;
  X >> c | b ;
]

routine = @quip [
  X >> 1             ;
  subroutine(1, 2, 3);
  subroutine(3, 2, 1);
]
```

<img src="https://github.com/npnguyen99/QViM/blob/main/assets/Subroutine.png?raw=true" height=300>

### 2.3. Measurements and classical bits
The virtual machine has a set of registers (can be customized) that stores classical bits (cBit). Access to these registers can be done using the syntax `C[index]`. cBits can be used as control bits to gate application. 
```julia
 @quip [ H >> 1 | (C[2], c[3]) ]
 ```

<img src="https://github.com/npnguyen99/QViM/blob/main/assets/cBit_Ctrl.png?raw=true" height=150>

 However, cBits and qBits can't be used together.
```julia
 @quip [ H >> 1 | (2, C[2]) ] # illegal, will error
```

A measurement will collapse a qBit to a definite state and write the result into a target cBit
```julia
@quip [ 1 => C[2] ]
```

<img src="https://github.com/npnguyen99/QViM/blob/main/assets/Measurement.png?raw=true" height=150>
