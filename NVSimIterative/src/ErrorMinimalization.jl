module ErrorMinimalization

export reconstructField

include("./NVSimulator.jl")
using .NVSimulator

using Unitful
using Measurements
using StaticArrays

# Numerical differentiation step size
# TODO: Make this user-controlled
const h = 1e-8
const hx = SA_F64[h, 0, 0]
const hy = SA_F64[0, h, 0]
const hz = SA_F64[0, 0, h]

# Coefficient for moving opposite to the gradient
# Higher values mean the algorithm converges faster, but setting it too high may cause pathological behaviour.
const step = 10

# Calculate error between measurement and simulation
function error(B, measured, Hnv)
    expected = expectedPeaks(Hnv, B)
    return sum(skipmissing((expected .- measured) .^ 2))
end

# Return RMSE of fit, taking into account missing values
function RMSE(B, measured, Hnv)
    err = error(B, measured, Hnv)u"GHz ^ 2"
    missingVals = count(ismissing, measured)
    return uconvert(u"MHz", sqrt(err/(length(measured) - missingVals)))
end

# Compute gradient numerically at a given point
function grad(func, x0)
    dx = (func(x0 + hx) - func(x0 - hx))/2h
    dy = (func(x0 + hy) - func(x0 - hy))/2h
    dz = (func(x0 + hz) - func(x0 - hz))/2h
    return SA_F64[dx, dy, dz]
end

# Field reconstruction algorithm
# measuredPeaks:        Sorted peak locations from measurement
# D:                    Calibrated D parameter (mean of the two zero field peaks)
# E:                    Calibrated E parameter (half-distance for the two zero field peaks)
# n:                    Number of iterations
# B0:                   Starting position for the search (vary this if you're getting nonsensical results)
function reconstructField(measuredPeaks; D = 2800.0, E = 0.0, n = 40_000, B0 = SA_F64[50,50,50])

    Hnv = calculateHnv(D, E)

    curriedError(B) = error(B, measuredPeaks, Hnv)
	
    Bs = fill(SA_F64[0,0,0], n)
    Bs[1] = B0
	
    errors = zeros(n)
    errors[1] = curriedError(B0)
	
    for i in 2:n
        gradient = grad(curriedError, Bs[i-1])
        Bs[i] = Bs[i-1] - step*gradient
        errors[i] = curriedError(Bs[i])
    end

    return (
        result = sort(Bs[end]), # The symmetries of the system make any permutation valid
        trajectory = Bs,
        errors = errors,
        RMSE = RMSE(Bs[end], measuredPeaks, Hnv)
    )
end

# Calculate the uncertainty of a single coordinate  
function coordinateUncertainty(peaks, i, u; D = 2800.0, E = 0.0, n = 40_000, B0 = SA_F64[50,50,50])
    function peakWrapper(p1, p2, p3, p4, p5, p6, p7, p8)
        pV = SA_F64[p1, p2, p3, p4, p5, p6, p7, p8]
        res = reconstructField(pV; D = D, E = E, n = n, B0 = B0).result
        return res[i]
    end
	
    local varied = peaks .± u
	
    # This is ugly, but the macro doesn't like vectors
    return @uncertain peakWrapper(
        varied[1], varied[2],
        varied[3], varied[4],
        varied[5], varied[6],
        varied[7], varied[8],
    ) 

end

function reconstructUncertainty(peaks, i, u; D = 2800.0, E = 0.0, n = 40_000, B0 = SA_F64[50,50,50])
	uB = [
		coordinateUncertainty(peaks, 1, u; D = D, E = E, n = n, B0 = B0),
		coordinateUncertainty(peaks, 2, u; D = D, E = E, n = n, B0 = B0),
		coordinateUncertainty(peaks, 3, u; D = D, E = E, n = n, B0 = B0),
	]u"Gauss"

	return uB
end

end
