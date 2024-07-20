module NVSim

using LinearAlgebra
using StaticArrays

include("Eigenstates.jl")
using .Eigenstates

orientation = normalize(SA_F64[0,0,-1])
Bs = range(0,800,1000)

res = eigenStates.(Ref(orientation), Bs)

end # module NVSim
