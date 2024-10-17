### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 0a0a94c8-653d-11ef-13fa-b98a27cb8524
using Plots

# ╔═╡ dbc2d3f9-5d70-4733-ae40-8750a8ea04ed
using CUDA

# ╔═╡ 77d5e493-08d0-4433-b8e5-0344de775b6e
using Peaks

# ╔═╡ c704e52b-e516-495c-84e1-3557466004b2
using JLD2

# ╔═╡ d67e3371-9b78-45dd-a2d0-38982080bca4
using StaticArrays

# ╔═╡ 7b7f582c-ce95-4f6a-a6de-29c3cf16e356
using LinearAlgebra

# ╔═╡ 0536ab28-19b0-4c68-ab5b-67f87cc5a2e2
using ProgressLogging

# ╔═╡ ad79aef6-907f-4c6d-8e15-e20ecf669e3b
begin
	const MHzToK = Float32(4.7991939324e-5)
	const GHzToK = Float32(4.7991939324e-2)
	const gToK = Float32(6.71714e-5) 
	const D = Float32(2870.2921719690626) # determined experimentally
	const E = Float32(2.8153124437499955) # determined experimentally
	const g = Float32(2.0026)
	md"Constants"
end

# ╔═╡ 6d1f7bba-54ea-47dd-a23b-e88185e8960b
begin
	const spinX = SMatrix{3,3,ComplexF32}([ 0 1 0
                                            1 0 1
                                            0 1 0 ] / sqrt(2))
	const spinY = SMatrix{3,3,ComplexF32}([  0 -im    0
                                             im   0  -im
                                              0   im   0 ] / sqrt(2))
	const spinZ = SMatrix{3,3,ComplexF32}([ 1  0  0
                                            0  0  0
                                            0  0 -1 ])
	md"Spin matrices"
end

# ╔═╡ 0c895b1e-feaa-4f3e-b6fe-fd4a4e054be4
begin
	const spinX2 = spinX^2
	const spinY2 = spinY^2
	const spinZ2 = spinZ^2
	md"Spin matrices squared"
end

# ╔═╡ 1cbe30b8-eb82-4bef-b774-0d19da33863a
begin
	const spinVector = @SVector [spinX, spinY, spinZ]
	md"Spin vector"
end

# ╔═╡ 27b097d2-9989-45d3-8725-3e9a4820024e
const Hnv = begin
    d = D * (spinZ2 - (spinX2 + spinY2 + spinZ2) / 3)
    e = E * (spinX2 - spinY2)
    MHzToK * (d + e)
end

# ╔═╡ e9b9d1bf-c257-45c5-93a2-25a99efaf73f
const dirs = normalize.([SA_F32[1,-1,-1], SA_F32[-1,1,-1], SA_F32[-1,-1,1], SA_F32[1,1,1]])

# ╔═╡ bc31145a-4d0c-41c9-ba8e-10af53d24bf5
md"#### Rotation functions"

# ╔═╡ f5859816-6659-4ce7-9e28-c20b803849ba
function Rx(t)
	SMatrix{3,3,Float32}(  1,       0,       0,
	                       0,  cos(t), -sin(t),
	                       0,  sin(t),  cos(t)  )
end

# ╔═╡ 2cddcd6c-f3d0-4308-9b60-b3c70507cab5
function Ry(t)
	SMatrix{3,3,Float32}(  cos(t),  0,  sin(t),
	                            0,  1,       0,
	                      -sin(t),  0,  cos(t)  )
end

# ╔═╡ 7a4581aa-14a0-467b-9ecf-1f351ce95e7b
function Rz(t)
	SMatrix{3,3,Float32}( cos(t), -sin(t),  0,
	                      sin(t),  cos(t),  0,
	                           0,       0,  1  )
end

# ╔═╡ 885a4f2a-fd99-478c-ab64-018660f620ce
function rotate(Rx, Ry, Rz, vect)
	Rz * Ry * Rx * vect
end

# ╔═╡ 327bf173-e909-4c0c-9cff-71efc2f6044e
function rotationMatrix(α, β, γ)
	Rz(γ) * Ry(β) * Rx(α)
end

# ╔═╡ f4733ba9-59a4-41d6-8690-dc3bcf9f7066
md"#### NV center simulation functions"

# ╔═╡ f86f903b-4418-4438-abca-afb8a8b02a9e
function unscaledZeeman(dir)
	gToK * g * sum(spinVector .* dir)
end

# ╔═╡ 0b624019-90c2-4790-b793-4d9d3f3833db
md"#### Input processing functions"

# ╔═╡ 9d37429c-a41b-4024-a02e-d934cd1ee48b
function findODMRPeaks(freqs, signal)
	pks, vals = findmaxima(signal)
	_, proms = peakproms!(pks, signal)
	decreasingProms = sort(proms, rev = true)
	peakFreqs = zeros(8)
	for i in 1:8
		peakIndex = findfirst(==(decreasingProms[i]), proms)
		peakFreqs[i] = freqs[pks[peakIndex]]
	end
	return sort(peakFreqs)
end

# ╔═╡ f8bf9242-71a0-4de7-b973-367dd1e42c45
begin
	testData = jldopen("../../data/2024-09-16.jld2")
	const peakVector = SVector{8, Float32}(testData["peaks"][4,:])
	VHs = testData["VHs"]
	IMs = testData["IMs"]
	md"Data import"
end

# ╔═╡ 05c0b368-be58-4ac5-ba7d-b917126b331a
# ╠═╡ disabled = true
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	peakMatrix = zeros(length(Bs), 8)
	for i in 1:length(Bs)
		peakMatrix[i,:] = findODMRPeaks(freqs, signals[i,:])
	end
	peakMatrix
end
  ╠═╡ =#

# ╔═╡ 23a93340-5f79-47ec-b135-0f77701306be
function calculateZeemans(Rx, Ry)
	rot = Rx * Ry
	return SVector(
		unscaledZeeman(rot * dirs[1]),
		unscaledZeeman(rot * dirs[2]),
		unscaledZeeman(rot * dirs[3]),
		unscaledZeeman(rot * dirs[4]),
	)
end

# ╔═╡ c2606b22-dbf1-4d58-8ca1-ea2e9f003615
function simKernel!(scores, Hzs, Bs)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
	j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
	k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

	if i <= size(scores)[1] && j <= size(scores)[2] && k <= size(scores)[3]
		@inbounds Hz = Hzs[i, j]
		@inbounds B = Bs[k]
		
		# Direction 1
		@inbounds vals1 = eigvals(Hnv + B * Hz[1])/GHzToK

		# Direction 2
		@inbounds vals2 = eigvals(Hnv + B * Hz[2])/GHzToK

		# Direction 3
		@inbounds vals3 = eigvals(Hnv + B * Hz[3])/GHzToK

		# Direction 4
		@inbounds vals4 = eigvals(Hnv + B * Hz[4])/GHzToK

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

		sortedPeaks = sort(peaks)

		for y in 1:8
			# Add the distance squared of the corresponding peaks to the score
			@inbounds scores[i, j, k] += (peakVector[y] - sortedPeaks[y])^2
		end
	end
    return
end

# ╔═╡ d44a7190-adb4-460a-949e-005d9afaed73
function fitODMR(αs, βs, Bs)
	n1 = length(αs)
	n2 = length(βs)
	n3 = length(Bs)
	
	dαs = CuArray{Float32}(αs)
	dβs = CuArray{Float32}(βs)

	Rxs = Rx.(dαs)
	Rys = Ry.(dβs)

	Hzs = calculateZeemans.(Rxs, Rys')

	dBs = CuArray{Float32}(Bs)
	
	threads = (8,8,4)
	blocks = ceil.(Int, (n1, n2, n3) ./ threads)

	scores = CUDA.zeros(Float32, n1, n2, n3)
	
	@cuda fastmath=true threads=threads blocks=blocks simKernel!(scores, Hzs, Bs)

	min = findmin(scores)
	bestfit = min[1]
	coords = min[2]

	α = αs[coords[1]]
	β = βs[coords[2]]
	B = Bs[coords[3]]
	
	return (
		bestfit = bestfit,
		orientation = (
			α = α,
			β = β,
			B = B
		)
	)

end

# ╔═╡ 35573e35-7bb0-4f06-9095-e2d0e4717ff6
# ╠═╡ disabled = true
#=╠═╡
fitODMR(range(0.00000001, π, 256), range(0.00000001, π, 256), range(200,215,256))
  ╠═╡ =#

# ╔═╡ 79dcf062-819b-4291-8f0a-b9355493261e
function segmentedFit(n, Bmin, Bmax)
	blockSize = 1014 # Batch size
	totalSize = blockSize * n

	bestFit = Inf
	bestOrientation = missing
	
	angles = range(0, π, totalSize + 1)[2:end]
	Bs = range(Bmin, Bmax, totalSize)

	@progress for i in 1:n, j in 1:n, k in 1:n
		iSlice = (((i - 1) * blockSize) + 1):(i * blockSize)
		jSlice = (((j - 1) * blockSize) + 1):(j * blockSize)
		kSlice = (((k - 1) * blockSize) + 1):(k * blockSize)

		fit = fitODMR(angles[iSlice], angles[jSlice], Bs[kSlice])

		if fit.bestfit < bestFit
			bestFit = fit.bestfit
			bestOrientation = fit.orientation
		end
	end

	return (bestfit = bestFit, orientation = bestOrientation)
		
end

# ╔═╡ 0857c507-69fb-4a13-9f9d-0aced7908381
CUDA.@profile segmentedFit(1, 200, 215)

# ╔═╡ 174177e9-f82d-4cf4-a396-bfea946d13b3
result = segmentedFit(2, 0, 500)

# ╔═╡ Cell order:
# ╠═0a0a94c8-653d-11ef-13fa-b98a27cb8524
# ╠═dbc2d3f9-5d70-4733-ae40-8750a8ea04ed
# ╠═77d5e493-08d0-4433-b8e5-0344de775b6e
# ╠═c704e52b-e516-495c-84e1-3557466004b2
# ╠═d67e3371-9b78-45dd-a2d0-38982080bca4
# ╠═7b7f582c-ce95-4f6a-a6de-29c3cf16e356
# ╠═0536ab28-19b0-4c68-ab5b-67f87cc5a2e2
# ╠═ad79aef6-907f-4c6d-8e15-e20ecf669e3b
# ╟─6d1f7bba-54ea-47dd-a23b-e88185e8960b
# ╟─0c895b1e-feaa-4f3e-b6fe-fd4a4e054be4
# ╟─1cbe30b8-eb82-4bef-b774-0d19da33863a
# ╟─27b097d2-9989-45d3-8725-3e9a4820024e
# ╟─e9b9d1bf-c257-45c5-93a2-25a99efaf73f
# ╟─bc31145a-4d0c-41c9-ba8e-10af53d24bf5
# ╟─f5859816-6659-4ce7-9e28-c20b803849ba
# ╟─2cddcd6c-f3d0-4308-9b60-b3c70507cab5
# ╟─7a4581aa-14a0-467b-9ecf-1f351ce95e7b
# ╟─885a4f2a-fd99-478c-ab64-018660f620ce
# ╟─327bf173-e909-4c0c-9cff-71efc2f6044e
# ╟─f4733ba9-59a4-41d6-8690-dc3bcf9f7066
# ╟─f86f903b-4418-4438-abca-afb8a8b02a9e
# ╟─0b624019-90c2-4790-b793-4d9d3f3833db
# ╟─9d37429c-a41b-4024-a02e-d934cd1ee48b
# ╠═f8bf9242-71a0-4de7-b973-367dd1e42c45
# ╠═05c0b368-be58-4ac5-ba7d-b917126b331a
# ╠═23a93340-5f79-47ec-b135-0f77701306be
# ╠═c2606b22-dbf1-4d58-8ca1-ea2e9f003615
# ╠═d44a7190-adb4-460a-949e-005d9afaed73
# ╠═35573e35-7bb0-4f06-9095-e2d0e4717ff6
# ╠═79dcf062-819b-4291-8f0a-b9355493261e
# ╠═0857c507-69fb-4a13-9f9d-0aced7908381
# ╠═174177e9-f82d-4cf4-a396-bfea946d13b3
