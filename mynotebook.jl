### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# ╔═╡ 49799ef0-ba51-11f0-28de-b37bf41e41e8
begin
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 8969a0b2-50f7-4573-9c63-40dcf7ef773e
using LinearAlgebra, Plots

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

# ╔═╡ 488991a5-2fda-4bc2-8702-74ad4459f1d4
## Constantes

begin
	μ=3.986004418e14 	# Paramètre gravitationnel terrestre
	Re=6371e3 			# Rayon moyen de la Terre
	ωe=7.2921150e-5 	# Vitesse de rotation de la Terre (considérée constante)

	nothing
end

# ╔═╡ e40e03da-6dad-40b0-b739-ae42d8ed98fb
## Structure pour stocker les satellites (définis par a, i, Ω et M0).

struct Sat
    a::Float64   # Demi-grand axe (m) → distance moyenne au centre de la Terre
    i::Float64   # Inclinaison orbitale (rad) → angle entre le plan orbital et l’équateur
    Ω::Float64   # Longitude du nœud ascendant (rad) → orientation du plan orbital autour de la Terre
    M0::Float64  # Anomalie moyenne initiale (rad) → position initiale du satellite sur son orbite
end


# ╔═╡ 2e3b8833-e011-41b1-be61-fdc9905e3115
## deg2rad et rad2deg functions

begin
	deg2rad(x)=x*pi/180
	rad2deg(x)=180*x/pi

	nothing
end

# ╔═╡ 03b47088-7ff3-4817-9012-7d156c638855
## Matrices de rotation

begin
	R1(θ)=[1 0 0; 0 cos(θ) -sin(θ); 0 sin(θ) cos(θ)] # Autour de l'axe X → Sert pour incliner le plan orbital d’un angle i (l'inclinaison).
	R3(θ)=[cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1] # Autour de l'axe Z → Sert pour tourner le plan orbital d’un angle Ω (ascension du nœud) ou de l’anomalie vraie ν.

	nothing
end

# ╔═╡ 326e41c9-3240-415a-a50a-135906a2e5c7
"""
Convertit un vecteur de position d’un satellite du repère ECI (Earth-Centered Inertial) vers le repère ECEF (Earth-Centered Earth-Fixed).
"""
function ecef_from_eci(r,t)
	return R3(-ωe*t)*r # page 129 - SE216 (v.Aout 2024)
end

# ╔═╡ d60e29d5-eeb4-4d19-87a9-56a0e62a8a75
"""
Convertit un vecteur position exprimé dans le repère ECEF (Earth-Centered Earth-Fixed) en latitude et longitude (degrés).
"""
function latlon_from_ecef(r)
    x, y, z = r # Vecteur position
    ρ = norm(r) # Distance entre le satellite et le centre de la Terre
	
    ϕ = asin(z / ρ) 	# Latitude -> angle entre le vecteur position et le plan équatorial.  
    λ = atan(y, x) 		# Longitude -> angle dans le plan équatorial entre l’axe X (Greenwich) et la projection du vecteur sur ce plan.

    return rad2deg(ϕ), rad2deg(λ)
end


# ╔═╡ 735be1be-0989-4e0d-8c36-910d80870534
"""
Calcule la position d’un satellite dans le repère inertiel (ECI) à un instant t, en supposant une orbite circulaire.
"""
function eci_pos(s::Sat, t)
    n = sqrt(μ / s.a^3) 			# Vitesse angulaire du satellite d'après la 3e loi de Kepler
    u = s.M0 + n * t 				# Position angulaire du satellite sur son orbite (Position initiale (M₀) + n * t)
    R = R3(s.Ω) * R1(s.i) * R3(u) 	# Matrice de transformation (Formule 7.35 - SE216 (v.Aout 2024))
    return R * [s.a, 0, 0] 			# Position finale du satellite
end


# ╔═╡ 48c2baeb-e62a-44a7-87d0-148705574470
"""
Génère une constellation de type Walker-Delta.\n
Utilisée pour répartir uniformément des satellites sur plusieurs plans orbitaux inclinés.
"""
function walker_delta(P, S, F, i_deg, a)
    sats = Sat[]
    for p in 0:P-1, s in 0:S-1 # On itère sur chaque plan orbital et sur chaque satellite par plan
        
		Ω = 2π * p / P # Ascension du nœud (Ω) → Chaque plan orbital est séparé de 360°/P autour de l’axe z.
		
        # Répartition uniforme des satellites sur chaque orbite (s/S),
        M0 = 2π * (s / S + (F * p) / (S * P)) # Anomalie moyenne initiale (M₀) + déphasage (F*p)/(S*P) pour éviter que les plans soient alignés.

        push!(sats, Sat(a, deg2rad(i_deg), Ω, M0)) # Création du satellite en fonction de a, i, Ω, M₀
    end
    return sats
end


# ╔═╡ 82695929-6414-4e51-ae34-c882211530c0
"""
Vérifie si un satellite est visible depuis un point au sol donné,
en tenant compte d’un angle d’élévation minimal (ε).
"""
function visible(r_ecef, lat_deg, lon_deg, eps_deg)
    ϕ = deg2rad(lat_deg)
    λ = deg2rad(lon_deg)

    g = Re * [cos(ϕ) * cos(λ), cos(ϕ) * sin(λ), sin(ϕ)] / Re # Vecteur qui pointe du centre de la Terre vers le point au sol défini par sa latitude ϕ et sa longitude λ. (repère ECEF)

    r = r_ecef 	# Vecteur position (ECEF)
    ρ = norm(r) # Distance entre le satellite et le centre de la Terre

    cosγ = dot(r, g) / (ρ * 1.0) 	# Produit scalaire normalisé
    γ = acos(clamp(cosγ, -1, 1))  	# Angle entre la direction du satellite et celle du point au sol
	
    ψmax = acos((Re / ρ) * cos(deg2rad(eps_deg))) # Angle au centre de la Terre jusqu’où le satellite reste visible au-dessus d’un certain angle d’élévation ε.

    return γ ≤ ψmax # Satellite visible avec au moins l'élévation minimale spécifiée
end


# ╔═╡ b04c08c9-62ad-4039-a90c-86f9e5478c55
"""
Calcule la fraction (%) de points visibles au moins par un satellite à un instant donné.
"""
function coverage_fraction(sats, t, latmin, latmax, eps_deg; dlat=2, dlon=2)
    pts = 0
    covered = 0

    r_ecef = [ecef_from_eci(eci_pos(s, t), t) for s in sats] # Position des satellites à l'instant t

    for lat in latmin:dlat:latmax, lon in -180:dlon:180
        pts += 1 
        for r in r_ecef 
            if visible(r, lat, lon, eps_deg)
                covered += 1 # +1 si un satellite voit le point 
                break
            end
        end
    end
    return covered / pts # Fraction couverte
end

# ╔═╡ 4899cb2f-04a9-40a6-9f56-0a63a8a37db3
"""
Affiche les zones couvertes par les satellites
"""
function show_coverage_heatmap!(sats, t, latmin, latmax, eps_deg; dlat=1, dlon=1)
    lats = collect(latmin:dlat:latmax)
    lons = collect(-180:dlon:180)
    r_ecef = [ecef_from_eci(eci_pos(s, t), t) for s in sats] # Position des satellites à l'instant t

    M = falses(length(lats), length(lons))

    for (i, lat) in enumerate(lats), (j, lon) in enumerate(lons)
        M[i, j] = any(r -> visible(r, lat, lon, eps_deg), r_ecef)
    end
    heatmap(lons, lats, Int.(M); aspect_ratio=1, colorbar=:false, ticks=:false, axis=:false)
end

# ╔═╡ 4d1173b5-3158-4535-abdd-14da09801f09
"""
Plot la Terre, les latitudes ±60°, et les positions instantanées des satellites en 3D.
"""
function plot_constellation!(sats, t; Rearth = Re)	
	# Hémisphère arrière de la Terre (longitudes ~ [-π, 0])
	u = range(-pi/2, pi/2, 60); v1 = range(-pi, 0, 60)
	xs = [Re*cos(ui)*cos(vi) for ui in u, vi in v1]
	ys = [Re*cos(ui)*sin(vi) for ui in u, vi in v1]
	zs = [Re*sin(ui)         for ui in u, vi in v1]
	p = surface(xs, ys, zs, color=:lightblue, opacity=0.03, linecolor=:transparent)
	
	# Hémisphère avant de la Terre (longitudes ~ [0, π])
	v2 = range(0, pi, 60)
	xs = [Re*cos(ui)*cos(vi) for ui in u, vi in v2]
	ys = [Re*cos(ui)*sin(vi) for ui in u, vi in v2]
	zs = [Re*sin(ui)         for ui in u, vi in v2]
	p = surface!(p, xs, ys, zs, color=:lightblue, opacity=0.05, linecolor=:transparent)

	θs = range(0, 2π, 100)
    for s in sats
        R = R3(s.Ω) * R1(s.i)  # orientation du plan orbital
        orb = [ecef_from_eci(R * [s.a * cos(θ), s.a * sin(θ), 0.0], t) for θ in θs]
        Xorb = [r[1] for r in orb]; Yorb = [r[2] for r in orb]; Zorb = [r[3] for r in orb]
    p = plot3d!(p, Xorb, Yorb, Zorb, lw=0.5, c=:gray, alpha=0.5)  # Tracé des plans orbitaux
    end
	
    r_ecef = [ecef_from_eci(eci_pos(s, t), t) for s in sats] # Position des satellites
    X = [r[1] for r in r_ecef]
    Y = [r[2] for r in r_ecef]
    Z = [r[3] for r in r_ecef]
    p = scatter3d!(p, X, Y, Z, marker=:circle, ms=3.5, color=:orange, aspect_ratio=:equal, legend=false, colorbar=false, size=(600,600)) # Plot des satellites

	return p
end

# ╔═╡ 6730fbd6-2cdd-4f88-a00c-182a601e6d97
## Paramètres

begin
	h_km=500.0 			# Altitude des satellites
	eps_deg=10 			# Elevation minimale
	i_deg=60.0 			# Inclinaison orbitale (angle avec l'équateur) pour savoir les lattitudes couvertes
	P=3 				# Nombre de plans orbitaux
	S=3 				# Nombre de satellites par plan
	F=0 				# Facteur de déphasage (pour décaler la position des satellites entre les plans afin d'éviter qu'ils soient alignés)
	a=Re + h_km*1e3 	# Demi-grand axe pour calculer la période orbitale
	t=0.0 				# Temps en secondes

	nothing
end

# ╔═╡ bfa07a7f-842b-4ea1-9cfb-ff902d694492
sats=walker_delta(P,S,F,i_deg,a)

# ╔═╡ 74249c3d-41d6-40dc-9f76-cc97d5df7ea0
cov = coverage_fraction(sats, t, -i_deg, i_deg, eps_deg; dlat=2, dlon=2)

# ╔═╡ 1332fe9d-cd3d-4c04-9009-42d2b80ddfb1
println("Couverture ±60° instantanée (ε=$(eps_deg)°) : $(round(cov*100,digits=2))% avec $(length(sats)) satellites.")

# ╔═╡ 27790501-c853-4b7d-9fdc-ab78585745fa
show_coverage_heatmap!(sats, t, -i_deg, i_deg, 0; dlat=1, dlon=1)

# ╔═╡ 171b6ed2-7180-4299-a6dc-12010c4c65c4
plot_constellation!(sats,t)

# ╔═╡ Cell order:
# ╟─ffc4fe26-29ab-405c-abb4-441945e251f0
# ╟─49799ef0-ba51-11f0-28de-b37bf41e41e8
# ╠═8969a0b2-50f7-4573-9c63-40dcf7ef773e
# ╠═488991a5-2fda-4bc2-8702-74ad4459f1d4
# ╠═e40e03da-6dad-40b0-b739-ae42d8ed98fb
# ╠═2e3b8833-e011-41b1-be61-fdc9905e3115
# ╠═03b47088-7ff3-4817-9012-7d156c638855
# ╟─326e41c9-3240-415a-a50a-135906a2e5c7
# ╟─d60e29d5-eeb4-4d19-87a9-56a0e62a8a75
# ╟─735be1be-0989-4e0d-8c36-910d80870534
# ╟─48c2baeb-e62a-44a7-87d0-148705574470
# ╟─82695929-6414-4e51-ae34-c882211530c0
# ╟─b04c08c9-62ad-4039-a90c-86f9e5478c55
# ╟─4899cb2f-04a9-40a6-9f56-0a63a8a37db3
# ╟─4d1173b5-3158-4535-abdd-14da09801f09
# ╠═6730fbd6-2cdd-4f88-a00c-182a601e6d97
# ╠═bfa07a7f-842b-4ea1-9cfb-ff902d694492
# ╠═74249c3d-41d6-40dc-9f76-cc97d5df7ea0
# ╠═1332fe9d-cd3d-4c04-9009-42d2b80ddfb1
# ╠═27790501-c853-4b7d-9fdc-ab78585745fa
# ╠═171b6ed2-7180-4299-a6dc-12010c4c65c4
