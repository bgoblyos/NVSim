### A Pluto.jl notebook ###
# v0.19.47

using Markdown
using InteractiveUtils

# ╔═╡ 157cec57-9f57-4bc1-9d6a-b593beb72dfb
using LinearAlgebra

# ╔═╡ 160880f2-95c9-11ef-2c25-6f273442942a
center = [622, 476]

# ╔═╡ 1e8077e8-b384-4d05-883b-1ff86f18a61c
before = [411, 312]

# ╔═╡ da448c97-993e-4dce-aba5-94c3631cb594
after = [272, 377]

# ╔═╡ a79209e9-816c-4e94-9729-74db06d22419
v1 = before - center

# ╔═╡ 632d5010-aede-4ae5-b35d-0da9bc23454b
v2 = after - center

# ╔═╡ 6e6330a7-5d68-4294-a698-880f59e67519
getAngle(a, b) = acos((a ⋅ b)/(norm(a)*norm(b)))

# ╔═╡ 7d013ca1-d960-4671-a652-5af4fdd21b4d
getAngle(v1, v2) / π * 180

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.1"
manifest_format = "2.0"
project_hash = "ac1187e548c6ab173ac57d4e72da1620216bce54"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

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

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"
"""

# ╔═╡ Cell order:
# ╠═157cec57-9f57-4bc1-9d6a-b593beb72dfb
# ╠═160880f2-95c9-11ef-2c25-6f273442942a
# ╠═1e8077e8-b384-4d05-883b-1ff86f18a61c
# ╠═da448c97-993e-4dce-aba5-94c3631cb594
# ╠═a79209e9-816c-4e94-9729-74db06d22419
# ╠═632d5010-aede-4ae5-b35d-0da9bc23454b
# ╠═6e6330a7-5d68-4294-a698-880f59e67519
# ╠═7d013ca1-d960-4671-a652-5af4fdd21b4d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
