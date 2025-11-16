### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 49799ef0-ba51-11f0-28de-b37bf41e41e8
begin
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 8969a0b2-50f7-4573-9c63-40dcf7ef773e
using LinearAlgebra, Plots

# ╔═╡ 6f253d7b-c6be-4755-a81e-75c8bd13c642
using Base.Threads # Pour la performance

# ╔═╡ daca6429-e148-42d3-9499-bfd1dbc04531
begin
	using Logging
	Logging.disable_logging(LogLevel(1000));
end

# ╔═╡ a123ac24-44c3-4c00-aed6-eed6cffe14f2
begin
	using Random
	Random.seed!(178);
end

# ╔═╡ bf43c4f4-11ce-455b-a3d2-5e1c11ab40d5
begin
    include("./SatsLEO/SatsLEO.jl")
    using .SatsLEO
end

# ╔═╡ 287f599c-0932-4d5c-ae28-16c6488d585a
using FileIO, PlutoUI

# ╔═╡ ffc4fe26-29ab-405c-abb4-441945e251f0
html"""
 <! -- this adapts the width of the cells to display its being used on -->
<style>
	main {
		margin: 0 auto;
		max-width: 2000px;
    	padding-left: max(160px, 10%);
    	padding-right: max(160px, 10%);
	}
</style>
"""

# ╔═╡ 6730fbd6-2cdd-4f88-a00c-182a601e6d97
## Paramètres

begin
	h_km=550 			# Altitude des satellites
	eps_deg=10 			# Elevation minimale nécéssaire pour voir le satellite depuis le sol
	i_deg=30 			# Inclinaison orbitale (angle avec l'équateur) pour savoir les lattitudes couvertes
	P=2 				# Nombre de plans orbitaux
	S=7 				# Nombre de satellites par plan
	F=0 				# Facteur de déphasage (pour décaler la position des satellites entre les plans afin d'éviter qu'ils soient alignés)
	a=Re + h_km*1e3 	# Demi-grand axe pour calculer la période orbitale
	t=0 				# Temps en secondes

	sats=walker_delta(P,S,F,i_deg,a)
end

# ╔═╡ 1ef59fe9-61f3-410a-b022-17d1592a2fc3
begin
	Ptest = 6
	Ntest = 21
	best_vec, best_cov = evolve_vec(Ptest, Ntest, F, i_deg, a, eps_deg; popsize=30, generations=40, Cmin=0.0)
	cov_final, _ = eval_constellation(best_vec, F, i_deg, a, eps_deg; n=100, dlat=1, dlon=1)
	best_vec, best_cov, cov_final
end

# ╔═╡ aa86b509-4bc9-4a95-84b6-9063f2e36b63
begin
    results = []
    Cmin = 95

    for P in 2:10, N in 2:25
        vec = random_vec(P, N)
        cov, Nt = eval_constellation(vec, F, i_deg, a, eps_deg)
        push!(results, (P=P, vec=vec, N=Nt, cov=cov))
    end

	good = filter(r -> r.cov ≥ Cmin, results)
	
	if isempty(good)
	    error("Aucune constellation ne dépasse Cmin = $Cmin")
	end
	
	Ns    = getfield.(good, :N)
	minN  = minimum(Ns)
	
	cands = filter(r -> r.N == minN, good)
	covs  = getfield.(cands, :cov)
	idx   = argmax(covs)
	
	best  = cands[idx]
	#good
	#cands
end

# ╔═╡ 41c3fdd3-e596-4f88-beda-16157552fea9
sats2 = myconstellation([8 8 8],1,i_deg,a)

# ╔═╡ 4985fe03-5e41-4ad6-899e-d864639107f8
coverage_fraction(sats, t, -i_deg, i_deg, eps_deg)

# ╔═╡ a3c17370-10bb-4f00-ad67-626b244318d6
coverage_fraction(sats2, t, -i_deg, i_deg, eps_deg)

# ╔═╡ 47fee197-c516-4624-a8a8-63345ade9841
begin
	p_1 = show_coverage_heatmap(sats2, t, eps_deg)
	p_2 = plot_constellation(sats2, t)
	plot(p_1, p_2; layout=(1,2), size=(1300,600))
end

# ╔═╡ f8f4c11b-b506-45df-b6b4-50abbe999c64
PlutoUI.LocalResource("./coverage_fraction_période.png")

# ╔═╡ 5b53534e-2155-4675-a271-454b92e3c96e
PlutoUI.LocalResource("./convergence_mean_coverage.png")

# ╔═╡ 84ff8c5f-5c92-4336-acb4-cd643d3e56e6
# Pour faire des GIFs sur CDN :)
if !isfile("./cov.gif") #!isfile("./sats.gif")
    folder = mktempdir()
	Tmax = 10000
	step = 100
    for t in 0:step:Tmax
		show_coverage_heatmap(sats, t, eps_deg)
		#plot_constellation(sats,t)
        savefig(joinpath(folder, "frame_$(Int(t/step)).png"))
    end
	
    frames = [load(joinpath(folder, "frame_$i.png")) for i in 0:Int(Tmax/step)]
    gr()
    save("./cov.gif", cat(frames..., dims=3))
	#save("./sats.gif", cat(frames_..., dims=3))
end

# ╔═╡ 89b95e5c-684c-44ca-9455-469e3bb97129
md"""
Click here to reload the GIF : $(@bind reload Button("Reload"))
"""

# ╔═╡ 8926723b-d835-4818-9073-89eea4b0dea4
begin
	reload
	PlutoUI.LocalResource("./sats.gif")
end

# ╔═╡ b2388c97-997d-4afd-a681-2b86b7c1458a
begin
	reload
	PlutoUI.LocalResource("./cov.gif")
end

# ╔═╡ Cell order:
# ╟─ffc4fe26-29ab-405c-abb4-441945e251f0
# ╟─49799ef0-ba51-11f0-28de-b37bf41e41e8
# ╠═8969a0b2-50f7-4573-9c63-40dcf7ef773e
# ╠═6f253d7b-c6be-4755-a81e-75c8bd13c642
# ╠═daca6429-e148-42d3-9499-bfd1dbc04531
# ╠═a123ac24-44c3-4c00-aed6-eed6cffe14f2
# ╠═bf43c4f4-11ce-455b-a3d2-5e1c11ab40d5
# ╠═1ef59fe9-61f3-410a-b022-17d1592a2fc3
# ╠═aa86b509-4bc9-4a95-84b6-9063f2e36b63
# ╠═6730fbd6-2cdd-4f88-a00c-182a601e6d97
# ╠═41c3fdd3-e596-4f88-beda-16157552fea9
# ╠═4985fe03-5e41-4ad6-899e-d864639107f8
# ╠═a3c17370-10bb-4f00-ad67-626b244318d6
# ╠═47fee197-c516-4624-a8a8-63345ade9841
# ╟─f8f4c11b-b506-45df-b6b4-50abbe999c64
# ╟─5b53534e-2155-4675-a271-454b92e3c96e
# ╠═287f599c-0932-4d5c-ae28-16c6488d585a
# ╠═84ff8c5f-5c92-4336-acb4-cd643d3e56e6
# ╟─89b95e5c-684c-44ca-9455-469e3bb97129
# ╠═8926723b-d835-4818-9073-89eea4b0dea4
# ╠═b2388c97-997d-4afd-a681-2b86b7c1458a
