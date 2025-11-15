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

# ╔═╡ 20df3e66-d9d0-43a9-9166-8c75cf70d8ac
# Bien pour initialiser un vecteur mais il faut l'améliorer
"""
balanced_vec(P, N)

Construit un vecteur de longueur P contenant une répartition aussi uniforme que possible
de N satellites entre les P plans orbitaux.
"""
function balanced_vec(P, N)
    base = N ÷ P              # Nombre minimal de satellites par plan
    r = N % P                 # Plans qui recevront un satellite supplémentaire
    v = fill(base, P)         # Répartition uniforme initiale
    for k in 1:r
        v[k] += 1             # Ajout des satellites restants
    end
    return v
end

# ╔═╡ 8345234a-9906-4182-9412-4c00dc6b703f
"""
improve_vec_fixed_N!(vec, F, i_deg, a, eps_deg)

Optimise la répartition de N satellites sur P plans orbitaux (nombre de satellites fixe) en déplaçant un satellite d’un plan vers un autre lorsque cela augmente la couverture moyenne sur une période. L’algorithme répète les déplacements les plus bénéfiques jusqu’à convergence ou jusqu’à 100 itérations.
"""
function improve_vec_fixed_N!(vec, F, i_deg, a, eps_deg)
	
    cov, N = eval_constellation(vec, F, i_deg, a, eps_deg)
    P = length(vec)

    for iter in 1:100
        best_gain = 0.0
        best_i, best_j = 0, 0

        for i in 1:P, j in 1:P
            i == j && continue
            vec[j] == 0 && continue

            vec[i] += 1; vec[j] -= 1 # Test de la config
	            cov2, _ = eval_constellation(vec, F, i_deg, a, eps_deg)
	            gain = cov2 - cov
	            if gain > best_gain # Si ça améliore, on garde i et j pour le refaire après
	                best_gain = gain
	                best_i, best_j = i, j
	            end
            vec[i] -= 1; vec[j] += 1 # On remet comme avant
        end

        if best_gain <= 0 # Si convergence, break
            break
        end

        vec[best_i] += 1; vec[best_j] -= 1 # On 
        cov, _ = eval_constellation(vec, F, i_deg, a, eps_deg)
    end

    vec, cov
end

# ╔═╡ 6730fbd6-2cdd-4f88-a00c-182a601e6d97
## Paramètres

begin
	h_km=500 			# Altitude des satellites
	eps_deg=10 			# Elevation minimale nécéssaire pour voir le satellite depuis le sol
	i_deg=30 			# Inclinaison orbitale (angle avec l'équateur) pour savoir les lattitudes couvertes
	P=2 				# Nombre de plans orbitaux
	S=7 				# Nombre de satellites par plan
	F=0 				# Facteur de déphasage (pour décaler la position des satellites entre les plans afin d'éviter qu'ils soient alignés)
	a=Re + h_km*1e3 	# Demi-grand axe pour calculer la période orbitale
	t=0 				# Temps en secondes

	sats=walker_delta(P,S,F,i_deg,a)
end

# ╔═╡ aa86b509-4bc9-4a95-84b6-9063f2e36b63
begin
    results = []
    Cmin = 0.95

    for P in 1:5, N in 2:21
        vec = balanced_vec(P, N)
        cov, Nt = eval_constellation(vec, F, i_deg, a, eps_deg)
        push!(results, (P=P, vec=vec, N=Nt, cov=cov))
    end

    good = filter(r -> r.cov ≥ Cmin, results)

    Ns   = getfield.(good, :N)
    idx  = argmin(Ns)
    best = good[idx]
end

# ╔═╡ b99bf2bc-6813-4a36-a55b-f9558ef2bab1
begin
	P_ = 8
	N_ = 18
	vec_ = balanced_vec(P_, N_)
	vec_opt_, cov_opt_ = improve_vec_fixed_N!(vec_, F, i_deg, a, eps_deg)
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
	p_2 = plot_constellation!(sats2, t)
	plot(p_1, p_2; layout=(1,2), size=(1300,600))
end

# ╔═╡ d97c7a1e-5929-4d50-aa6e-c1b029be1861
# ╠═╡ disabled = true
#=╠═╡
let
    h_km = 800
    eps_deg = 10
    i_deg = 30
    P = 2
    S = 7
    F = 0
    a = Re + h_km*1e3

	sats=walker_delta(P,S,F,i_deg,a)
	
	Ts = 0:100:100000
	covs = [coverage_fraction(sats, T, -90, 90, eps_deg)*100 for T in Ts]
	
	plot(Ts, round.(covs; digits=2),
	    xlabel = "Temps (s)",
	    ylabel = "Couverture (%)",
	    legend = false,
	    title = "Évolution de la couverture en fonction du temps",
	    markersize = 3,
	    color = :blue)
end
  ╠═╡ =#

# ╔═╡ 5b53534e-2155-4675-a271-454b92e3c96e
# ╠═╡ disabled = true
#=╠═╡
let
    h_km = 800
    eps_deg = 10
    i_deg = 30
    P = 2
    S = 5
    F = 0
    a = Re + h_km*1e3

    sats = walker_delta(P, S, F, i_deg, a)

    N = 200
    ns = 2 .* (1:N)
    vals = zeros(N)

    for k in 1:N
        vals[k] = mean_coverage_fraction(sats, -i_deg, i_deg, eps_deg; n=ns[k])
    end

    plot(
        ns, vals;
        xlabel = "Nombre d'échantillons temporels",
        ylabel = "Couverture moyenne (%)",
        title = "Convergence de la couverture moyenne\nWalker-Delta P=$(P), S=$(S), i=$(i_deg)°",
        lw = 2,
        markershape = :circle,
        markerstrokewidth = 0,
        legend = false,
        grid = true,
    )
end
  ╠═╡ =#

# ╔═╡ e67f0e51-c074-4ed6-b38a-fd4afbdaa9be
begin
	sats_ = myconstellation([1,1,10],1,i_deg,a)
	r_ecef = [[[ecef_from_eci(eci_pos(s, t),t) for s in sats_],coverage_fraction(sats_, t, -i_deg, i_deg, eps_deg)] for t in 1:10]
	#show_coverage_heatmap(sats_,t,eps_deg)
end

# ╔═╡ 27790501-c853-4b7d-9fdc-ab78585745fa
begin
	p1 = show_coverage_heatmap(sats, t, eps_deg)
	p2 = plot_constellation!(sats, t)
	plot(p1, p2; layout=(1,2), size=(1300,600))
end

# ╔═╡ 84ff8c5f-5c92-4336-acb4-cd643d3e56e6
# Pour faire des GIFs sur CDN :)
if !isfile("./cov.gif") #!isfile("./sats.gif")
    folder = mktempdir()
	Tmax = 10000
	step = 100
    for t in 0:step:Tmax
		show_coverage_heatmap(sats, t, eps_deg)
		#plot_constellation!(sats,t)
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
	#PlutoUI.LocalResource("./sats.gif")
end

# ╔═╡ b2388c97-997d-4afd-a681-2b86b7c1458a
begin
	reload
	#PlutoUI.LocalResource("./cov.gif")
end

# ╔═╡ Cell order:
# ╟─ffc4fe26-29ab-405c-abb4-441945e251f0
# ╟─49799ef0-ba51-11f0-28de-b37bf41e41e8
# ╠═8969a0b2-50f7-4573-9c63-40dcf7ef773e
# ╠═6f253d7b-c6be-4755-a81e-75c8bd13c642
# ╠═daca6429-e148-42d3-9499-bfd1dbc04531
# ╠═bf43c4f4-11ce-455b-a3d2-5e1c11ab40d5
# ╟─20df3e66-d9d0-43a9-9166-8c75cf70d8ac
# ╠═8345234a-9906-4182-9412-4c00dc6b703f
# ╠═aa86b509-4bc9-4a95-84b6-9063f2e36b63
# ╠═b99bf2bc-6813-4a36-a55b-f9558ef2bab1
# ╠═6730fbd6-2cdd-4f88-a00c-182a601e6d97
# ╠═41c3fdd3-e596-4f88-beda-16157552fea9
# ╠═4985fe03-5e41-4ad6-899e-d864639107f8
# ╠═a3c17370-10bb-4f00-ad67-626b244318d6
# ╠═47fee197-c516-4624-a8a8-63345ade9841
# ╟─d97c7a1e-5929-4d50-aa6e-c1b029be1861
# ╟─5b53534e-2155-4675-a271-454b92e3c96e
# ╠═e67f0e51-c074-4ed6-b38a-fd4afbdaa9be
# ╟─27790501-c853-4b7d-9fdc-ab78585745fa
# ╠═287f599c-0932-4d5c-ae28-16c6488d585a
# ╠═84ff8c5f-5c92-4336-acb4-cd643d3e56e6
# ╟─89b95e5c-684c-44ca-9455-469e3bb97129
# ╠═8926723b-d835-4818-9073-89eea4b0dea4
# ╠═b2388c97-997d-4afd-a681-2b86b7c1458a
