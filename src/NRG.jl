module NRG

using LinearAlgebra
using SparseArrays

#export plusTwo
#plusTwo(x) = return x+2

include("WilsonParam.jl")
include("FractalFillings.jl")
include("NRGiterations.jl")
include("Observables.jl")

end 
