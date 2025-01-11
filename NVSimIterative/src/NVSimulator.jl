module NVSimulator

export calculateHnv
export expectedPeaks

using StaticArrays
using LinearAlgebra

# Constants
const MHzToK = Float64(4.7991939324e-5)
const GHzToK = Float64(4.7991939324e-2)
const gToK = Float64(6.71714e-5) 
const g = Float64(2.0026)

# Spin matrices
const spinX = SMatrix{3,3,ComplexF64}([ 0 1 0
                                        1 0 1
                                        0 1 0 ] / sqrt(2))
const spinY = SMatrix{3,3,ComplexF64}([  0 -im    0
                                        im   0  -im
                                         0  im    0 ] / sqrt(2))
const spinZ = SMatrix{3,3,ComplexF64}([ 1  0  0
                                        0  0  0
                                        0  0 -1 ])

# Spin matrices squared
const spinX2 = spinX^2
const spinY2 = spinY^2
const spinZ2 = spinZ^2

# Spin vector
const spinVector = @SVector [spinX, spinY, spinZ]

# NV directions
const dirs = normalize.([SA_F64[1,-1,-1], SA_F64[-1,1,-1], SA_F64[-1,-1,1], SA_F64[1,1,1]])

# Calculate the ZFS Hamiltonian based on the supplied D and E parameters
function calculateHnv(D, E)
    d = D * (spinZ2 - (spinX2 + spinY2 + spinZ2) / 3)
    e = E * (spinX2 - spinY2)
    MHzToK * (d + e)
end

# Rotation matrix to bring v1 to v2
function R(v1, v2)
    u = v1 × v2
    c = v1 ⋅ v2
    s = norm(u)
    u = normalize(u)

    asym = SMatrix{3,3,Float32}(
            0, -u[3],  u[2],
         u[3],     0, -u[1],
        -u[2],  u[1],     0,
    )
	
    rotation_matrix = c * I + s * asym + (1 - c) * (u * u')
end

# Calculate Zeeman splitting Hamiltonian for unit field strength
function unscaledZeeman(dir)
    gToK * g * sum(spinVector .* dir)
end

# Compute all 4 Hamiltonians corresponding to the NV directions for a given magnetic field
function hamiltonians(Hnv, Bx, By, Bz)
    BVec = SA_F64[Bx, Bz, By]
    B = norm(BVec)
    rot = R(BVec/B, SA_F64[0, 0, 1])
    return SVector(
        unscaledZeeman(rot * dirs[1]),
        unscaledZeeman(rot * dirs[2]),
        unscaledZeeman(rot * dirs[3]),
        unscaledZeeman(rot * dirs[4]),
    ) * B .+ Ref(Hnv)
end

# Calculate the expected peaks corresponding to an induction vector
function expectedPeaks(Hnv, B)
	
    Hs = hamiltonians(Hnv, B...)
		
    # Direction 1
    @inbounds vals1 = eigvals(Hs[1])/GHzToK

    # Direction 2
    @inbounds vals2 = eigvals(Hs[2])/GHzToK

    # Direction 3
    @inbounds vals3 = eigvals(Hs[3])/GHzToK

    # Direction 4
    @inbounds vals4 = eigvals(Hs[4])/GHzToK

    @inbounds peaks = SVector(
        abs(vals1[1] - vals1[3]),
        abs(vals1[1] - vals1[2]),
        abs(vals2[1] - vals2[3]),
        abs(vals2[1] - vals2[2]),
        abs(vals3[1] - vals3[3]),
        abs(vals3[1] - vals3[2]),
        abs(vals4[1] - vals4[3]),
        abs(vals4[1] - vals4[2])
    )

    return sort(peaks)
end

end
