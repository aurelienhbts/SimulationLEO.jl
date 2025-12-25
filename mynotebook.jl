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
using LinearAlgebra, Plots, FileIO, PlutoUI

# ╔═╡ 6f253d7b-c6be-4755-a81e-75c8bd13c642
using Base.Threads # Pour la performance

# ╔═╡ daca6429-e148-42d3-9499-bfd1dbc04531
begin
	using Logging
	Logging.disable_logging(LogLevel(1)) # Enlever les warnings
end

# ╔═╡ bf43c4f4-11ce-455b-a3d2-5e1c11ab40d5
begin
    include("./SatsLEO/SatsLEO.jl") # Module SatsLEO dev pour le projet
    using .SatsLEO
end

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

# ╔═╡ 89b95e5c-684c-44ca-9455-469e3bb97129
md"""
Click here to reload the GIF : $(@bind reload Button("Reload"))
"""

# ╔═╡ b2388c97-997d-4afd-a681-2b86b7c1458a
begin
	reload
	PlutoUI.LocalResource("./Figures/cov.gif")
end

# ╔═╡ 8926723b-d835-4818-9073-89eea4b0dea4
begin
	reload
	PlutoUI.LocalResource("./Figures/sats.gif")
end

# ╔═╡ 7a7ad281-b620-4fcf-9daa-2d4fc2985a80
begin
	reload
	PlutoUI.LocalResource("./Figures/satstest.gif")
end

# ╔═╡ 58ae7adf-6552-434e-a086-8db4cf360c40
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

# ╔═╡ 84ff8c5f-5c92-4336-acb4-cd643d3e56e6
# Pour faire des GIFs sur CDN :)
if !isfile("./Figures/cov.gif") #!isfile("./Figures/sats.gif")
    folder = mktempdir()
	vec = [7 0 4 4 4 0]
	satsgif = myconstellation(vec, F, i_deg, a)
	Tmax = 10000
	step = 100
    Threads.@threads for t in 0:step:Tmax
		show_coverage_heatmap(satsgif, t, eps_deg)
		#plot_constellation(satsgif,t)
        savefig(joinpath(folder, "frame_$(Int(t/step)).png"))
    end
	
    frames = [load(joinpath(folder, "frame_$i.png")) for i in 0:Int(Tmax/step)]
    gr()
    save("./Figures/cov.gif", cat(frames..., dims=3))
	#save("./Figures/sats.gif", cat(frames..., dims=3))
end

# ╔═╡ 1003088c-a8da-4f5b-a00c-e86adc808559
# Pour faire des GIFs sur CDN :)
if !isfile("./Figures/satstest.gif")
    folder_ = mktempdir()
	vec_ = [5 5 5 5 5 5 5 5 5 5]
	satsgif_ = myconstellation(vec_, F, 40, a)
	Tmax_ = 10000
	step_ = 100
    Threads.@threads for t_ in 0:step_:Tmax_
		plot_constellation(satsgif_,t_)
        savefig(joinpath(folder_, "frame_$(Int(t_/step_)).png"))
    end
	
    frames_ = [load(joinpath(folder_, "frame_$i.png")) for i in 0:Int(Tmax_/step_)]
    gr()
    save("./Figures/satstest.gif", cat(frames_..., dims=3))
end

# ╔═╡ 47fee197-c516-4624-a8a8-63345ade9841
# ╠═╡ disabled = true
#=╠═╡
begin
	sats2 = walker_delta(2, 10, 1, i_deg ,a)
	p_1 = show_coverage_heatmap(sats2, t, eps_deg)
	p_2 = plot_constellation(sats2, t)
	plot(p_1, p_2; layout=(1,2), size=(1300,600))
end
  ╠═╡ =#

# ╔═╡ f8f4c11b-b506-45df-b6b4-50abbe999c64
PlutoUI.LocalResource("./Figures/coverage_fraction_période.png")

# ╔═╡ 5b53534e-2155-4675-a271-454b92e3c96e
PlutoUI.LocalResource("./Figures/convergence_mean_coverage.png")

# ╔═╡ c08f4857-5666-4346-b66d-7ed4a943973d
empty!(FITCACHE)

# ╔═╡ fff4d4f2-3965-4dbc-a76b-e649827a063d
let
	iter = 5
	a_test = 500 *1e3 + Re
	i_deg_test = 63.4
	Pmax = 6
	Ninit = 20
	configs = Dict()
	for _ in 1:iter
		best_vec, cov_final, N_final = evolve_vec(Pmax, Ninit, F, i_deg_test, a_test, eps_deg; nsats=1, generations=100, Ctarget=95, Ncoef=0.0, K=0.0)
		if !haskey(configs,(best_vec, cov_final, N_final))
			configs[(best_vec, cov_final, N_final)] = 1/iter * 100
		else 
			configs[(best_vec, cov_final, N_final)] += 1/iter * 100
		end
	end
	configs = sort!(collect(configs), by = x -> x[1][2], rev=true)
end

# ╔═╡ 2f6f8114-02b6-4380-8b4c-dc77a7862880
# p_add et p_rem variables en fonction de Ctarget et best_cov (ok)

# ╔═╡ 05899639-28f6-43e4-81ab-7f68e72ada48
# Ajouter un Nmax (ok)

# ╔═╡ 2bbfeca6-9f80-46a6-a2cd-e6c2a8fb4bf1
# Ajouter le fait qu'il faille plusieurs satellites par location (ok)
# Verifier 24 pour galileo 
# Voir l'impact si il y en a qui est defect
# Comparer avec des constellations actuelles

# ╔═╡ b269560d-ec58-4a02-a50a-5f3a6f91a111
# PDOP → qualité du signal gnss (dilution of precision)
# Calculer la line of sight des satellites 
# Si la valeur est tres élevée, le gps n'est pas très précis

# ╔═╡ f9e6e837-1a81-41b4-8c6a-5bbfcb44c122
# Equations très connues concernant la communication
# Si l'amplitude du signal niveau donné → constel optimale

# ╔═╡ 402783a9-e443-40ee-a1ba-93c7899d9085
# Comment gérer le fait qu'il soit possible que plusieurs vec soient ajoutés au FITCACHE en même temps ?

# ╔═╡ 89fbe235-53ba-4c9b-83a9-2fbc8fea3396
PlutoUI.LocalResource("./concurrent_write_dict.png")

# ╔═╡ a38a8820-06d1-4716-af41-48bddf0af9b8
# copie dans chaque thread et puis recombiner
# Plus facile à dire qu'à faire

# ╔═╡ 7779dcf1-0f43-4eb5-92b8-c3b039c3f760
empty!(FITCACHE_fixedN)

# ╔═╡ 6c8b3b55-bfc2-43d7-9d88-f6119ff61ee1
let
	iter = 10
	atest = 550 *1e3 + Re
	Pmax = 6
	Ntest = 20
	configs = Dict()
	for _ in 1:iter
		best_vec, best_cov = evolve_vec_fixedN(Pmax, Ntest, F, i_deg, atest, eps_deg; popsize=20, generations=50, Cmin=90, Pbonus=true)
		if !haskey(configs,(best_vec,best_cov))
			configs[(best_vec,best_cov)] = 1/iter * 100
		else 
			configs[(best_vec,best_cov)] += 1/iter * 100
		end
	end
	configs = sort!(collect(configs), by = x -> x[1][2], rev=true)
end

# ╔═╡ 556db4d2-7329-4c7d-b1ee-d98245d8756a
PlutoUI.LocalResource("./Figures/convergence_cov_altitude.png")

# ╔═╡ e6121f95-8f68-47e3-97ae-0c6422b80ee4
PlutoUI.LocalResource("./Figures/convergence_cov_altitude_2.png")

# ╔═╡ Cell order:
# ╟─ffc4fe26-29ab-405c-abb4-441945e251f0
# ╟─49799ef0-ba51-11f0-28de-b37bf41e41e8
# ╠═8969a0b2-50f7-4573-9c63-40dcf7ef773e
# ╠═6f253d7b-c6be-4755-a81e-75c8bd13c642
# ╠═daca6429-e148-42d3-9499-bfd1dbc04531
# ╠═bf43c4f4-11ce-455b-a3d2-5e1c11ab40d5
# ╠═84ff8c5f-5c92-4336-acb4-cd643d3e56e6
# ╟─89b95e5c-684c-44ca-9455-469e3bb97129
# ╟─b2388c97-997d-4afd-a681-2b86b7c1458a
# ╟─8926723b-d835-4818-9073-89eea4b0dea4
# ╠═1003088c-a8da-4f5b-a00c-e86adc808559
# ╟─7a7ad281-b620-4fcf-9daa-2d4fc2985a80
# ╠═58ae7adf-6552-434e-a086-8db4cf360c40
# ╠═47fee197-c516-4624-a8a8-63345ade9841
# ╟─f8f4c11b-b506-45df-b6b4-50abbe999c64
# ╟─5b53534e-2155-4675-a271-454b92e3c96e
# ╠═c08f4857-5666-4346-b66d-7ed4a943973d
# ╠═fff4d4f2-3965-4dbc-a76b-e649827a063d
# ╠═2f6f8114-02b6-4380-8b4c-dc77a7862880
# ╠═05899639-28f6-43e4-81ab-7f68e72ada48
# ╠═2bbfeca6-9f80-46a6-a2cd-e6c2a8fb4bf1
# ╠═b269560d-ec58-4a02-a50a-5f3a6f91a111
# ╠═f9e6e837-1a81-41b4-8c6a-5bbfcb44c122
# ╠═402783a9-e443-40ee-a1ba-93c7899d9085
# ╟─89fbe235-53ba-4c9b-83a9-2fbc8fea3396
# ╠═a38a8820-06d1-4716-af41-48bddf0af9b8
# ╠═7779dcf1-0f43-4eb5-92b8-c3b039c3f760
# ╠═6c8b3b55-bfc2-43d7-9d88-f6119ff61ee1
# ╟─556db4d2-7329-4c7d-b1ee-d98245d8756a
# ╟─e6121f95-8f68-47e3-97ae-0c6422b80ee4
