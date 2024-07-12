### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ dc663401-3356-4cd0-8e01-7619a0e7df42
begin
using PlutoUI
import PlutoUI:combine
import PlutoUI:ExperimentalLayout
end

# ╔═╡ 0ff62a0d-53e6-48b9-8246-3f7b89499417
using TTVFaster,Test

# ╔═╡ ebea622c-b4a2-4fdd-bb3e-2925695da863
using PyPlot, DelimitedFiles

# ╔═╡ acc61da3-1bf6-4413-a3a9-3bacfd318def
@bind reset_to_default_values Button("Reset Sliders")

# ╔═╡ d08ce186-d8e5-4515-9b46-9bb62260936c
begin
function conditions(system::String,iplanet,values::Vector{T}) where (T<:Real)
	# Create sliders for given system, planet id (number or letter)
	# Here, masses must be given in log
	# 		this allows more precise selection of slider value
	# MSUN= 1.98892e33; MEARTH = 5.9742e27 #cgs
	names=["log_mu","per","trans0","ecosw","esinw"]
	min_vals=[log(1e-8),-2*values[2],-500,-1,-1]
	max_vals=[log(1e-2), 2*values[2], 1.05*values[3], 1, 1]
	return combine() do Child
		inputs=[[
			md"""
			 $(names[i]): $(Child(names[i], Slider(min_vals[i]:0.01:max_vals[i],default=values[i])))	
			"""
			for i=1:3
		];[md"""
			 $(names[i]): $(Child(names[i], Slider(min_vals[i]:0.001:max_vals[i],default=values[i])))	
			""" 
		for i=4:5]]
		md"""
		###### $(system) $(iplanet)
		$(inputs)
			"""
	end
	return values
end
# Allow for input of planet_plane structure from TTVFaster
function conditions(system_plname::String,p::Planet_plane_hk{T}) where (T<:Real)
	names=["mu","period","trans0","ecosw","esinw"]
	min_vals=[1, 1, -200, -1, -1]
	max_vals=[50, 500,5000, 1, 1]
	return combine() do Child
		inputs=[[
			md"""
			 $("mu"): $(Child("mu", Slider(min_vals[1]:0.01:max_vals[1],default=p.mass_ratio)))	
			"""];[md"""
			 $("period"): $(Child("period", Slider(min_vals[2]:0.001:max_vals[2],default=p.period)))	
			"""];[md"""
			 $("trans0"): $(Child("trans0", Slider(min_vals[3]:0.001:max_vals[3],default=p.trans0)))	
			"""];[md"""
			 $("ecosw"): $(Child("ecosw", Slider(min_vals[4]:0.001:max_vals[4],default=p.ecosw)))	
			"""];[md"""
			 $("esinw"): $(Child("esinw", Slider(min_vals[4]:0.001:max_vals[4],default=p.esinw)))	
			""" ]]
		md"""
		###### $(system_plname)
		$(inputs)
			"""
	end
	return values
end
md""" Create sliders based on orbital elements of a planet
"""
end

# ╔═╡ 4702121d-b071-4dd3-a43b-00fdf9c33f70
begin
data=readdlm("../examples/trappist1efg_planets.txt",',',Float64)
# Define initial conditions (i.e. the default parameters) for sliders
# Bind sliders to parameters that we want to model, add bonds as necessary
system="Trappist-1"
# For demonstration purposes only,
# this does not include the effects from other planets in the system
p1=@bind planet1 conditions(system,"e",
				[log(data[1]) ;data[2:5]])
p2=@bind planet2 conditions(system,"f",
				[log(data[6]) ;data[7:10]])
p3=@bind planet3 conditions(system,"g",
				[log(data[11]) ;data[12:15]])
# Bind reset button to sliders
let reset_to_default_values
	p1
	p2
	p3
end
# Create grid layout for multiple sliders, add bonds as necessary
# Layout must be in matrix form (i.e. no empty elements)
function create_grid(bonds)
	bond1,bond2,bond3=bonds[1:end]
	grid=PlutoUI.ExperimentalLayout.grid([
		bond1 bond2 bond3 
		])
	return grid
end
md" Approximate values for $system planets e/f/g"
end

# ╔═╡ 3c6c85c7-299a-4d5f-8769-a236e31bfd14
create_grid([p1,p2,p3])

# ╔═╡ b16d8dd8-00c3-4bc0-8a93-bedd356dc577
param_sliders=[
	planet1.log_mu planet1.per planet1.trans0 planet1.ecosw planet1.esinw
	planet2.log_mu planet2.per planet2.trans0 planet2.ecosw planet2.esinw
	planet3.log_mu planet3.per planet3.trans0 planet3.ecosw planet3.esinw
] # Binds parameter sliders that we want to plot

# ╔═╡ 6ec3cebb-f763-431e-be22-2f7ce3072fc0
# View individual parameters for a planet
planet1

# ╔═╡ 7cb10c9e-24eb-48b8-9b35-697fec3ec79c
begin 
	# from test_ttv.jl
	include("test_ttv.jl")
	hk1=TTVFaster.Planet_plane_hk(data[1], data[2], data[3], data[4], data[5])
	hk2=TTVFaster.Planet_plane_hk(data[6], data[7], data[8], data[9], data[10])
function kepler62_test()
    data=readdlm("../examples/kepler62ef_planets.txt",',',Float64)
    @time ttv1,ttv2=test_ttv(5,40,20,data,WriteOutput=false,num_evals=100000);
    inner=readdlm("../examples/inner_ttv.txt")
    outer=readdlm("../examples/outer_ttv.txt")
	inner_ref=readdlm("inner_ttv.txt.ref")
	outer_ref=readdlm("outer_ttv.txt.ref")
	diffs=[maximum(abs.(inner_ref .- inner)) maximum(abs.(outer_ref .- outer))]
	max_diff=maximum(diffs)
    return round(max_diff)
end
	
@testset "TTVFaster.jl" begin
    @test kepler62_test() == (0.0)
	# @test trappist1_test()
end
end

# ╔═╡ 525e9b6c-8255-42d7-a704-3924f0d27213
begin
# Plot the provided system parameters, given transit counts
# params must be an matrix with shape (nplanet, 5)
function plot_TTVs(ntrans,params)
	# @assert(length(ntrans)==nplanet)
	@assert(size(params)[1]==length(ntrans))
	nplanet=size(params)[1]
	
	# Reshape parameters into array with shape(1 x 5*nplanet)
	new_params=reshape(transpose(params),1,nplanet*5)
	# Get pair-wise TTVs
	function decompose_ttvs(nplanet,ntrans,params)
		jmax = 5
		pair_ttvs = zeros(nplanet,nplanet,maximum(ntrans))
		for i=1:nplanet-1,j=i+1:nplanet
			param = [exp(params[(i-1)*5+1]);params[(i-1)*5+2:i*5] ;exp(params[(j-1)*5+1]); params[(j-1)*5+2:j*5]]
			ttv = TTVFaster.ttv_nplanet(2,jmax,[ntrans[i];ntrans[j]],param)
			pair_ttvs[i,j,1:ntrans[i]] = ttv[1,1:ntrans[i]] #planet i wrt planet j
			pair_ttvs[j,i,1:ntrans[j]] = ttv[2,1:ntrans[j]] #planet j wrt planet i
		end
		return pair_ttvs
	end
	pair_ttvs=decompose_ttvs(nplanet,ntrans,new_params) .* (24 * 60)
	# Assuming that there's 2 transiting planets
	n1=ntrans[1];n2=ntrans[2]
	t1  = collect(new_params[3] .+ new_params[2] .* range(0,stop=n1-1,length=n1)) 
	t2  = collect(new_params[7] .+ new_params[6] .* range(0,stop=n2-1,length=n2))
	p1_ttvs=sum([pair_ttvs[1,iplanet,1:n1] for iplanet=1:nplanet],dims=1)[1]
	p2_ttvs=sum([pair_ttvs[2,iplanet,1:n2] for iplanet=1:nplanet],dims=1)[1]
	# Plot TTVs
	fig=plt.figure(figsize=(6,4))
	ax=fig.add_subplot(221)
	ax3=fig.add_subplot(222,sharey=ax,sharex=ax)
	ax3.plot(t1,transpose(p1_ttvs),color="grey")
	ax3.set_title("Total TTVs")
	ax2=fig.add_subplot(223)
	ax4=fig.add_subplot(224,sharey=ax2,sharex=ax2)
	ax4.plot(t2,transpose(p2_ttvs),color="grey")
	for iplanet=1:nplanet 
		ax2.plot(t2,pair_ttvs[2,iplanet,1:n2],label=string("planet ",iplanet))
		ax.plot(t1,pair_ttvs[1,iplanet,1:n1])#,label=string("planet",iplanet))
	end
	ax2.set_ylabel("TTVs [min]")
	ax.set_ylabel("TTVs [min]")
	ax2.set_xlabel(string("Time [JD ",L"$- t_{0,1}$","]"))
	ax4.set_xlabel(string("TIme [JD ",L"$- t_{0,1}$","]"))
	ax.minorticks_on() ; ax2.minorticks_on() ; plt.tight_layout()
	fig.legend(loc="upper left",title="Source of Perturbations",fontsize="small",ncol=nplanet)
	return fig
end
end

# ╔═╡ 4db7da85-ae23-4a34-bb69-5aab907d7145
plot_TTVs([40,20,2],param_sliders)

# ╔═╡ e00727ea-f5fb-49bc-83c2-21fb50dbe7e1
# using HypertextLiteral:@htl

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PyPlot = "d330b81b-6aea-500a-939a-2ce795aea3ee"
TTVFaster = "d84f081e-b698-44a3-a477-911041168508"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[compat]
PlutoUI = "~0.7.59"
PyPlot = "~2.11.5"
TTVFaster = "~0.2.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "b19db3927f0db4151cb86d073689f2428e524576"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.10.2"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

[[PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "9816a3826b0ebf49ab4926e2b18842ad8b5c8f04"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.96.4"

[[PyPlot]]
deps = ["Colors", "LaTeXStrings", "PyCall", "Sockets", "Test", "VersionParsing"]
git-tree-sha1 = "0371ca706e3f295481cbf94c8c36692b072285c2"
uuid = "d330b81b-6aea-500a-939a-2ce795aea3ee"
version = "2.11.5"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TTVFaster]]
deps = ["DelimitedFiles"]
git-tree-sha1 = "7fea3b01e0d667fd77281ad215efeff82948bbf3"
uuid = "d84f081e-b698-44a3-a477-911041168508"
version = "0.2.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═acc61da3-1bf6-4413-a3a9-3bacfd318def
# ╠═3c6c85c7-299a-4d5f-8769-a236e31bfd14
# ╠═4db7da85-ae23-4a34-bb69-5aab907d7145
# ╠═b16d8dd8-00c3-4bc0-8a93-bedd356dc577
# ╠═4702121d-b071-4dd3-a43b-00fdf9c33f70
# ╠═6ec3cebb-f763-431e-be22-2f7ce3072fc0
# ╠═7cb10c9e-24eb-48b8-9b35-697fec3ec79c
# ╠═d08ce186-d8e5-4515-9b46-9bb62260936c
# ╠═525e9b6c-8255-42d7-a704-3924f0d27213
# ╠═dc663401-3356-4cd0-8e01-7619a0e7df42
# ╠═e00727ea-f5fb-49bc-83c2-21fb50dbe7e1
# ╠═0ff62a0d-53e6-48b9-8246-3f7b89499417
# ╠═ebea622c-b4a2-4fdd-bb3e-2925695da863
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
