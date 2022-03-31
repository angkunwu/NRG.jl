module NRG

using LinearAlgebra

#greet() = print("Hello World!")

export plusTwo

plusTwo(x) = return x+2

include("WilsonParam.jl")


end # module
