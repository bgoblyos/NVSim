module Eigenstates

using LinearAlgebra
using StaticArrays

export eigenStates

const MHzToK = 4.7991939324e-5
const GHzToK = 4.7991939324e-2
const gToK = 6.71714e-5

const D = 2875
const E = 0
const g = 2.0026

const spinX = SMatrix{3,3}([ 0 1 0
                             1 0 1
                             0 1 0 ] / sqrt(2))

const spinY = SMatrix{3,3}([  0 -im   0
                             im   0 -im
                              0   im   0 ] / sqrt(2))

const spinZ = SMatrix{3,3}([ 1  0  0
                             0  0  0
                             0  0 -1 ])

spinVector = @SVector [spinX, spinY, spinZ]

spinX2 = spinX^2
spinY2 = spinY^2
spinZ2 = spinZ^2

# ZFS Hamiltonian
Hnv = begin
    d = D * (spinZ2 - (spinX2 + spinY2 + spinZ2)/3)
    e = E * (spinX2 - spinY2)
    MHzToK * (d + e)
end

function eigenStates(dir, B)
    Hzeeman = B * gToK * g * sum(spinVector .* dir)
    H = Hzeeman + Hnv
    eigenvals = eigvals(H)
    return eigenvals / GHzToK
end


end
