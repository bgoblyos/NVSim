### A Pluto.jl notebook ###
# v0.19.47

using Markdown
using InteractiveUtils

# ╔═╡ ef95f9ba-95c3-11ef-185b-193fa9811b10
using LinearAlgebra

# ╔═╡ 6ddb2a92-b65e-417e-b3ff-7c192707b766
using StaticArrays

# ╔═╡ 8e49ac22-a705-4c1d-a63f-1fedd2809223
using Combinatorics

# ╔═╡ de808a72-2b70-43db-8add-12710ae58f71
# ╠═╡ disabled = true
#=╠═╡
function C3(axis, θ)
	# using Rodrigues' rotation formula
	x = axis[1]
	y = axis[2]
	z = axis[3]
	
	K = SMatrix{3,3}([  0 -z  y
                        z  0 -x
                       -y  x  0 ])

	R = I + sin(θ)*K + (1 - cos(θ)) * K^2
	
end
  ╠═╡ =#

# ╔═╡ 8f850873-7018-4fee-bc1b-7f4b2a6fa074
const dirs = normalize.([SA_F32[1,-1,-1], SA_F32[-1,1,-1], SA_F32[-1,-1,1], SA_F32[1,1,1]])

# ╔═╡ a1df879e-b611-4f3a-9b13-db1a80e205dd
function shift(v)
	return [v[3], v[1], v[2]]
end

# ╔═╡ 609839c3-ef24-4386-b99c-4fa4fc230817
getAngle(a, b) = acos((a ⋅ b)/(norm(a)*norm(b)))

# ╔═╡ de322617-50b5-4430-a356-cbf6c815606f
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

# ╔═╡ 92c1ed7f-9d18-4638-8f68-fd829c4bc76e
function allAngles(B1, B2)
	B1s = allCombinations(B1)
	unique(sort(getAngle.(Ref(B2), B1s)))
end

# ╔═╡ d0866a6d-8a29-416a-a92c-822699212588
md"""
## Experimental data
"""

# ╔═╡ 6bb4c4df-df30-4247-8ecd-ffd784e71e66
B1 = [0.828665, 0.440783, 0.345000]

# ╔═╡ fe24c7bc-6a66-43c7-bed5-c466d385a17f
B2 = [0.803962, 0.206103, 0.557823 ]

# ╔═╡ e796914a-719f-4822-887e-6cd79f0ced7a
allAngles(B1, B2) / π * 180

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Combinatorics = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[compat]
Combinatorics = "~1.0.2"
StaticArrays = "~1.9.7"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.1"
manifest_format = "2.0"
project_hash = "af61c377fee52b80a3cfbb108c4575130d49211b"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
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

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

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

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"
"""

# ╔═╡ Cell order:
# ╠═ef95f9ba-95c3-11ef-185b-193fa9811b10
# ╠═6ddb2a92-b65e-417e-b3ff-7c192707b766
# ╠═8e49ac22-a705-4c1d-a63f-1fedd2809223
# ╠═de808a72-2b70-43db-8add-12710ae58f71
# ╠═8f850873-7018-4fee-bc1b-7f4b2a6fa074
# ╠═a1df879e-b611-4f3a-9b13-db1a80e205dd
# ╠═609839c3-ef24-4386-b99c-4fa4fc230817
# ╠═de322617-50b5-4430-a356-cbf6c815606f
# ╠═92c1ed7f-9d18-4638-8f68-fd829c4bc76e
# ╟─d0866a6d-8a29-416a-a92c-822699212588
# ╠═6bb4c4df-df30-4247-8ecd-ffd784e71e66
# ╠═fe24c7bc-6a66-43c7-bed5-c466d385a17f
# ╠═e796914a-719f-4822-887e-6cd79f0ced7a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
