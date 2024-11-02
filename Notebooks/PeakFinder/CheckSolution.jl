### A Pluto.jl notebook ###
# v0.19.47

using Markdown
using InteractiveUtils

# ╔═╡ 831832ba-961a-11ef-37ab-c1da01352df0
using StaticArrays

# ╔═╡ b25f24a2-502b-46e2-8311-3205803a967e
using LinearAlgebra

# ╔═╡ b98178cf-39fa-405a-8fea-c4b91cb593ba
using JLD2

# ╔═╡ fa4e6cfa-40cc-4cae-b9df-67dd9ae9d574
using Combinatorics

# ╔═╡ 11247173-0499-4208-a04c-89b4f3266650
begin
	testData = jldopen("../../data/2024-09-09.jld2")
	const peakVector = SVector{8, Float64}(testData["peaks"][2,:])
	md"Data import"
end

# ╔═╡ 50b4d218-3223-4e63-a39d-99a31b063c8f
begin
	const MHzToK = Float64(4.7991939324e-5)
	const GHzToK = Float64(4.7991939324e-2)
	const gToK = Float64(6.71714e-5) 
	const D = Float64(2870.292)
	const E = Float64(2.815312)
	const g = Float64(2.0026)
	md"Constants"
end

# ╔═╡ 4e9bddb3-1f63-4213-b599-3d80e86de57b
begin
	const spinX = SMatrix{3,3,ComplexF64}([ 0 1 0
                                            1 0 1
                                            0 1 0 ] / sqrt(2))
	const spinY = SMatrix{3,3,ComplexF64}([  0 -im    0
                                             im   0  -im
                                              0   im   0 ] / sqrt(2))
	const spinZ = SMatrix{3,3,ComplexF64}([ 1  0  0
                                            0  0  0
                                            0  0 -1 ])
	md"Spin matrices"
end

# ╔═╡ 8f231eb7-1f3b-4957-b7e5-03a0fee8b93a
begin
	const spinX2 = spinX^2
	const spinY2 = spinY^2
	const spinZ2 = spinZ^2
	md"Spin matrices squared"
end

# ╔═╡ 46039e49-47b6-4bf3-80a1-fbeaf0065890
begin
	const spinVector = @SVector [spinX, spinY, spinZ]
	md"Spin vector"
end

# ╔═╡ 24793d04-5cc2-4da3-8d94-ec08272fbf4b
const Hnv = begin
    d = D * (spinZ2 - (spinX2 + spinY2 + spinZ2) / 3)
    e = E * (spinX2 - spinY2)
    MHzToK * (d + e)
end

# ╔═╡ f8b90923-081f-40c4-83ba-1f6ee7486e39
const dirs = normalize.([SA_F64[1,-1,-1], SA_F64[-1,1,-1], SA_F64[-1,-1,1], SA_F64[1,1,1]])

# ╔═╡ 8db341ef-1a81-4fa3-a8af-e61445bd59a5
function unscaledZeeman(dir)
	gToK * g * sum(spinVector .* dir)
end

# ╔═╡ ae6d7813-f5a3-40f1-86f6-ec89c6efb6fc
function R(v1, v2) # rotation matrix that brings v1 to v2
	u = v1 × v2
	c = v1 ⋅ v2
	s = norm(u)
	u = normalize(u)

	asym = [
             0    -u[3]  u[2]
             u[3]  0    -u[1]
            -u[2]  u[1]  0
        ]
	
    rotation_matrix = c * I + s * asym + (1 - c) * (u * u')
end

# ╔═╡ cb037454-0dc7-485f-8358-9298158fe584
function calculateZeemans(θ, ϕ)
	Bvec = SA_F32[
		sin(θ) * cos(ϕ),
		sin(θ) * sin(ϕ),
		cos(θ)
	]
	rot = R(u, SA_F32[0, 0, 1])
	return SVector(
		unscaledZeeman(rot * dirs[1]),
		unscaledZeeman(rot * dirs[2]),
		unscaledZeeman(rot * dirs[3]),
		unscaledZeeman(rot * dirs[4]),
	)
end

# ╔═╡ 97ecbe27-0876-4258-bdb5-b5fea6b0cda8
function cpuSim(Bvec)
	B = norm(Bvec)
	u = normalize(Bvec)

	Hz = calculateZeemans(u)

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

	return sum((peakVector .- sortedPeaks).^2)
end

# ╔═╡ 3eff6168-e901-48bf-9538-40f9b9a68271
function shift(v)
	return [v[3], v[1], v[2]]
end

# ╔═╡ bf87c94d-987e-4138-a7f8-c1b51c188bda
function allCombinations(B0)
	Bperm = permutations(B0)
	results = []
	for B in Bperm
		for i in [1, -1], j in [1, -1], k in [1, -1]
			D = Diagonal([i,j,k])
			push!(results, D*B)
		end
	end	
	results
end

# ╔═╡ 0a0e4364-12ca-4681-bcb9-0f1e228eb352
B0 = 102.368*[0.803962, 0.206103, 0.557823 ]

# ╔═╡ 69c10eef-bdac-414e-938a-44e0682f30b3
Bs = allCombinations(B0)

# ╔═╡ 3fdfb564-8890-4a06-a8ad-06dd19306b9b
cpuSim.(Bs)

# ╔═╡ 68f0ce1e-c546-4335-9a0d-66a07ba71a6f
R([1,0,0], [0,0,1]) * [1,0,0]

# ╔═╡ 14847ffe-8fc0-42e7-b13d-f70e54ea4df2
Bi = -102.368*[0.557823, 0.206103, 0.803962]

# ╔═╡ 8366b08d-565d-4f27-a42e-a6a260b3644a
cpuSim(Bi)

# ╔═╡ 651855ba-2bdd-4ad9-9498-84c2652eff53
unique(allCombinations([1,2,3]))

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Combinatorics = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[compat]
Combinatorics = "~1.0.2"
JLD2 = "~0.5.5"
StaticArrays = "~1.9.7"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.1"
manifest_format = "2.0"
project_hash = "8a3f500a73600eb6f40a65ec6f4336cef21f26b2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "62ca0547a14c57e98154423419d8a342dca75ca9"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.4"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "PrecompileTools", "Requires", "TranscodingStreams"]
git-tree-sha1 = "aeab5c68eb2cf326619bf71235d8f4561c62fe22"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.5.5"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "eeafab08ae20c62c44c8399ccb9354a04b80db50"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.7"

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

    [deps.StaticArrays.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╠═831832ba-961a-11ef-37ab-c1da01352df0
# ╠═b25f24a2-502b-46e2-8311-3205803a967e
# ╠═b98178cf-39fa-405a-8fea-c4b91cb593ba
# ╠═fa4e6cfa-40cc-4cae-b9df-67dd9ae9d574
# ╠═11247173-0499-4208-a04c-89b4f3266650
# ╠═50b4d218-3223-4e63-a39d-99a31b063c8f
# ╠═4e9bddb3-1f63-4213-b599-3d80e86de57b
# ╠═8f231eb7-1f3b-4957-b7e5-03a0fee8b93a
# ╠═46039e49-47b6-4bf3-80a1-fbeaf0065890
# ╠═24793d04-5cc2-4da3-8d94-ec08272fbf4b
# ╠═f8b90923-081f-40c4-83ba-1f6ee7486e39
# ╠═8db341ef-1a81-4fa3-a8af-e61445bd59a5
# ╠═ae6d7813-f5a3-40f1-86f6-ec89c6efb6fc
# ╠═cb037454-0dc7-485f-8358-9298158fe584
# ╠═97ecbe27-0876-4258-bdb5-b5fea6b0cda8
# ╠═3eff6168-e901-48bf-9538-40f9b9a68271
# ╠═bf87c94d-987e-4138-a7f8-c1b51c188bda
# ╠═0a0e4364-12ca-4681-bcb9-0f1e228eb352
# ╠═69c10eef-bdac-414e-938a-44e0682f30b3
# ╠═3fdfb564-8890-4a06-a8ad-06dd19306b9b
# ╠═68f0ce1e-c546-4335-9a0d-66a07ba71a6f
# ╠═14847ffe-8fc0-42e7-b13d-f70e54ea4df2
# ╠═8366b08d-565d-4f27-a42e-a6a260b3644a
# ╠═651855ba-2bdd-4ad9-9498-84c2652eff53
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
