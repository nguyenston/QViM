# QViM: Quantum Virtual Machine
As the name suggests this package allows for the creation of simulated quantum machines.

## 1. Basic Usage

A quantum circuit can be thought of as a series of gate applications.
Similarly, the virtual machine class `QVM` is essentially a turing machine with an internal state (wave function) that takes in
a list of `QInstr` (QInstr Program or QuIP) and execute on it procedually.

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
The expression below (reads **apply `U` on `a` given `b`**) represents the following circuit:
```julia
@quip [U >> a | b]
```

<img src="https://github.com/npnguyen99/QViM/blob/main/assets/Basic_QInstr.png?raw=true" height=150>

In the case of multiple target or control qBits, use comma separated tuples.
```julia
@quip [F >> (a1, a2) | (b1, b2)]
```

<img src="https://github.com/npnguyen99/QViM/blob/main/assets/Multibit_QInstr.png?raw=true" height=350>

### 2.2. Control on zero
Putting the `!` operator in front of a control bit to indicate control on zero on that bit
```julia
@quip [ H >> 1 | !2 ]
```

<img src="https://github.com/npnguyen99/QViM/blob/main/assets/Ctrl_on_zero.png?raw=true" height=150>

### 2.3. Function call
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

### 2.4. Measurements and classical bits
The virtual machine has a set of registers (can be customized) that stores classical bits (cBit). Access to these registers can be done using the syntax `C[index]`. cBits can be used as and act like control bits in gate application. 
```julia
 @quip [ H >> 1 | (C[2], C[3]) ]
 ```

<img src="https://github.com/npnguyen99/QViM/blob/main/assets/cBit_Ctrl.png?raw=true" height=150>

 However, cBits and qBits can't be used together.
```julia
 @quip [ H >> 1 | (2, C[2]) ] # illegal, will error
```

A measurement will collapse a qBit to a definite state and write the result into a target cBit.
```julia
@quip [ 1 => C[2] ]
```

<img src="https://github.com/npnguyen99/QViM/blob/main/assets/Measurement.png?raw=true" height=150>

### 2.5. List comprehension, conditional expression, and nesting QuIPs
Conveniently, the `@quip` macro can parse list comprehension expressions. Bellow is an example of its usage in combination with a conditional expression.
```julia
@quip [i == 3 ? 
          H >> i : 
          X >> i | (i + 1) 
       for i in 1:3]
```
**Note:** When using non-atomic expressions like `i + 1` for qBits, surround them with parantheses `(i + 1)`.

<img src="https://github.com/npnguyen99/QViM/blob/main/assets/List_Comp_example.png?raw=true" height=300>

Since the `@quip` macro flattens the array expression, these example below are valid.
```julia
@quip [
  X >> 1 ;
  [
    H >> 3 ;
    X >> 3 ;
  ];
  Y >> 2 ;
]

@quip [
  H >> 1 ;
  [X >> (i + 1) | i for i in 1:3] ;
  3 => C[1] ;
]
```
## 3. The QVM object

### 3.1. Constructor
The constructor of the `QVM` object takes the form of:
```julia
QVM(quip::Vector{QInstr}, nbit::Int;  # required parameters
            dict=GATES, nreg=32)      # optional keyword parameters
```

Where:
* `quip` is the main program to be run
* `nbit` is the number of qBits simulated
* `dict` is the set of primitive operations that can be applied on the wave function
* `nreg` is the number of classical registers

### 3.2. Relevant functions
* `showstate(vm::QVM, n=2; io::IO=stdout)`: pretty print the current state of `vm`, like stack frames and current instruction
* `step!(vm::QVM)`: step through a single instruction
* `execute!(vm::QVM; verbose::Bool=false, io::IO=stdout)`: step through all instructions to the end of the main program
