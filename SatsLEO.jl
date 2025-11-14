module SatsLEO

using LinearAlgebra, Plots

export  μ, Re, ωe, Sat, deg2rad, rad2deg, R1, R3, 
        ecef_from_eci, latlon_from_ecef, eci_pos, walker_delta, myconstellation,
        visible, coverage_fraction, mean_coverage_fraction, show_coverage_heatmap,
        plot_constellation!, plot_earth

## Constantes
const μ = 3.986004418e14 	# Paramètre gravitationnel terrestre
const Re = 6371e3 		    # Rayon moyen de la Terre
const ωe = 7.2921150e-5 	# Vitesse de rotation de la Terre (considérée constante)

## Structure pour stocker les satellites (définis par a, i, Ω et M0).
struct Sat
    a::Float64   # Demi-grand axe (m) → distance moyenne au centre de la Terre
    i::Float64   # Inclinaison orbitale (rad) → angle entre le plan orbital et l’équateur
    Ω::Float64   # Longitude du nœud ascendant (rad) → orientation du plan orbital autour de la Terre
    M0::Float64  # Anomalie moyenne initiale (rad) → position du satellite sur son orbite
end

## deg2rad et rad2deg functions
deg2rad(x)=x*pi/180
rad2deg(x)=180*x/pi

## Matrices de rotation
R1(θ)=[1 0 0; 0 cos(θ) -sin(θ); 0 sin(θ) cos(θ)] # Autour de l'axe X → Sert pour incliner le plan orbital d’un angle i (l'inclinaison).
R3(θ)=[cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1] # Autour de l'axe Z → Sert pour tourner le plan orbital d’un angle Ω (ascension du nœud) ou de l’anomalie vraie ν.

"""
ecef_from_eci(r,t)

Convertit un vecteur de position d’un satellite du repère ECI (Earth-Centered Inertial) vers le repère ECEF (Earth-Centered Earth-Fixed).
"""
function ecef_from_eci(r,t)
	return R3(-ωe*t)*r # page 129 - SE216 (v.Aout 2024)
end

"""
latlon_from_ecef(r)

Convertit un vecteur position exprimé dans le repère ECEF (Earth-Centered Earth-Fixed) en latitude et longitude (degrés).
"""
function latlon_from_ecef(r)
    x, y, z = r # Vecteur position
    ρ = norm(r) # Distance entre le satellite et le centre de la Terre
	
    ϕ = asin(z / ρ) 	# Latitude -> angle entre le vecteur position et le plan équatorial.  
    λ = atan(y, x) 		# Longitude -> angle dans le plan équatorial entre l’axe X (Greenwich) et la projection du vecteur sur ce plan.

    return rad2deg(ϕ), rad2deg(λ)
end

"""
eci_pos(sat,t)

Calcule la position d’un satellite dans le repère inertiel (ECI) à un instant t, en supposant une orbite circulaire.
"""
function eci_pos(sat,t;mu=μ)
    n = sqrt(mu / sat.a^3) 			    # Vitesse angulaire du satellite d'après la 3e loi de Kepler
    u = sat.M0 + n * t 				    # Position angulaire du satellite sur son orbite (Position initiale (M₀) + n * t)
    R = R3(sat.Ω) * R1(sat.i) * R3(u) 	# Matrice de transformation (Formule 7.35 - SE216 (v.Aout 2024))
    return R * [sat.a, 0, 0] 			# Position finale du satellite
end

"""
walker_delta(P,S,F,i_deg,a)

Initialise une constellation de type Walker-Delta. Utilisée pour répartir uniformément des satellites sur plusieurs plans orbitaux inclinés.
"""
function walker_delta(P,S,F,i_deg,a)
    sats = Sat[]
    for p in 0:P-1, s in 0:S-1 # On itère sur chaque plan orbital et sur chaque satellite par plan
        
		Ω = 2π * p / P # Ascension du nœud (Ω) → Chaque plan orbital est séparé de 360°/P autour de l’axe z.
		
        # Répartition uniforme des satellites sur chaque orbite (s/S),
        M0 = 2π * (s / S + (F * p) / (S * P)) # Anomalie moyenne initiale (M₀) + déphasage (F*p)/(S*P) pour éviter que les plans soient alignés.

        push!(sats, Sat(a, deg2rad(i_deg), Ω, M0)) # Création du satellite en fonction de a, i, Ω, M₀
    end
    return sats
end

"""
myconstellation(vec,i_deg,a)

Initialise une constellation de type Walker-delta mais avec plus de flexibilité. vec est un vecteur comportant le nombre de satellite par plan orbital.
"""
function myconstellation(vec,F,i_deg,a)
	sats = Sat[]
	P = length(vec)
	for p in 1:P
		S = vec[p]
		for s in 1:S
			Ω = 2π * p / P # Ascension du nœud (Ω) → Chaque plan orbital est séparé de 360°/P autour de l’axe z.
			
			# Répartition uniforme des satellites sur chaque orbite (s/S),
        	M0 = 2π * (s / S + (F * p) / (S * P)) # Anomalie moyenne initiale (M₀) + déphasage (F*p)/(S*P) pour éviter que les plans soient alignés.
			
			push!(sats, Sat(a, deg2rad(i_deg), Ω, M0)) # Création du satellite en fonction de a, i, Ω, M₀
		end
	end
	return sats
end

"""
visible(r_ecef, lat_deg, lon_deg, eps_deg)

Vérifie si un satellite est visible depuis un point au sol donné,
en tenant compte d’un angle d’élévation minimal (ε).
"""
function visible(r_ecef, lat_deg, lon_deg, eps_deg)
    ϕ = deg2rad(lat_deg)
    λ = deg2rad(lon_deg)
	
	# Vecteur qui pointe du centre de la Terre vers le point au sol défini par sa latitude et sa longitude. (repère ECEF)
    gx = cos(ϕ) * cos(λ)
    gy = cos(ϕ) * sin(λ)
    gz = sin(ϕ)

    ρ = norm(r_ecef) # Distance entre le satellite et le centre de la Terre

    cosγ = (r_ecef[1]*gx + r_ecef[2]*gy + r_ecef[3]*gz) / ρ
    cosψmax = (Re / ρ) * cos(deg2rad(eps_deg)) # Angle entre la direction du satellite et celle du point au sol

    return cosγ >= cosψmax
end

"""
coverage_fraction(sats, t, latmin, latmax, eps_deg)

Calcule la fraction (%) de points visibles au moins par un satellite **à un instant t donné**.
Tient compte de l'élévation minimale nécessaire pour voir les satellites depuis le sol.
"""
function coverage_fraction(sats, t, latmin, latmax, eps_deg)
    r_ecef = [ecef_from_eci(eci_pos(s, t), t) for s in sats] # Position des satellites à l'instant t
    lats = collect(latmin:2:latmax)
    lons = collect(-180:2:180)

    nthreads = Threads.nthreads()
    covered_threads = fill(0, nthreads)
    pts_threads = fill(0, nthreads)

    Threads.@threads for i in eachindex(lats)
        tid = Threads.threadid()
        lat = lats[i]
        @inbounds for lon in lons
            pts_threads[tid] += 1
            for r in r_ecef
                if visible(r, lat, lon, eps_deg)
                    covered_threads[tid] += 1
                    break
                end
            end
        end
    end

    covered = sum(covered_threads)
    pts = sum(pts_threads)
    return 100 * covered / pts
end

"""
mean_coverage_fraction(sats, latmin, latmax, eps_deg)

Calcule la fraction (%) des points visibles au moins par un satellite en moyenne sur une période orbitale.
Tient compte de l'élévation minimale nécessaire pour voir les satellites depuis le sol (ε).
"""
function mean_coverage_fraction(sats, latmin, latmax, eps_deg; n=100)
    a = sats[1].a
    T = 2π * sqrt(a^3 / μ) # Periode orbitale en utilisant la 3e loi de Kepler
    ts = range(0, T; length=n) # Pour évaluer la couverture à 100 différents moments
    s = 0.0
    for t in ts
        s += coverage_fraction(sats, t, latmin, latmax, eps_deg)
    end
    return s / length(ts)
end

"""
show_coverage_heatmap(sats,t,eps_deg)

Affiche les zones couvertes par les satellites. Tient compte de l'élévation minimale nécessaire pour voir les satellites depuis le sol.
"""
function show_coverage_heatmap(sats,t,eps_deg)
    lats = collect(-90:1:90)
    lons = collect(-180:1:180)
    r_ecef = [ecef_from_eci(eci_pos(s, t), t) for s in sats]

    M = falses(length(lats), length(lons))

    for (i, lat) in enumerate(lats), (j, lon) in enumerate(lons)
        M[i, j] = any(r -> visible(r, lat, lon, eps_deg), r_ecef)
    end
	p = heatmap(lons, lats, Int.(M))
	p = plot!(p, xlabel = "Longitude [°]", ylabel = "Latitude [°]", title = "Zones couvertes par $(length(sats)) satellites")
	
    return plot!(p, aspect_ratio=1, colorbar=false, framestyle=:none)
end

"""
plot_constellation!(sats,t;Rearth=Re)

Plot les positions instantanées des satellites en 3D autour de la Terre.
"""
function plot_constellation!(sats,t;Rearth=Re)
	p = plot_earth()
	
	θs = range(0, 2π, 100)
    for s in sats
        R = R3(s.Ω) * R1(s.i)  # orientation du plan orbital
        orb = [R * [s.a * cos(θ), s.a * sin(θ), 0.0] for θ in θs]
        Xorb = [r[1] for r in orb]; Yorb = [r[2] for r in orb]; Zorb = [r[3] for r in orb]
    p = plot3d!(p, Xorb, Yorb, Zorb, lw=0.5, c=:gray, alpha=0.5)  # Tracé des plans orbitaux
    end
	
    r_ecef = [eci_pos(s, t) for s in sats] # Position des satellites
    X = [r[1] for r in r_ecef]
    Y = [r[2] for r in r_ecef]
    Z = [r[3] for r in r_ecef]
    p = scatter3d!(p, X, Y, Z, marker=:circle, ms=3.5, color=:orange, title="Visualisation des satellites en t = $(Int(t))s") # Plot des satellites

	return p
end

"""
plot_earth(;Rearth=Re)

Plot une sphère (la Terre) avec des lignes pour les lattitudes et l'équateur en bleu.
"""
function plot_earth(;Rearth=Re)
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

	# Ajout de lignes pour les lattitudes
	for lat_deg in -80:10:-10
        lat_rad = deg2rad(lat_deg)
        r = Rearth * cos(lat_rad)
        z = Rearth * sin(lat_rad)
        θs = range(0, 2π, 100)
        x = [r * cos(θ) for θ in θs]
        y = [r * sin(θ) for θ in θs]
        z_line = fill(z, length(θs))
        p = plot3d!(p, x, y, z_line, lw=0.3, c=:black, alpha=0.3)
    end
	for lat_deg = 0
        lat_rad = deg2rad(lat_deg)
        r = Rearth * cos(lat_rad)
        z = Rearth * sin(lat_rad)
        θs = range(0, 2π, 100)
        x = [r * cos(θ) for θ in θs]
        y = [r * sin(θ) for θ in θs]
        z_line = fill(z, length(θs))
        p = plot3d!(p, x, y, z_line, lw=0.3, c=:blue)
    end
	for lat_deg in 10:10:80
        lat_rad = deg2rad(lat_deg)
        r = Rearth * cos(lat_rad)
        z = Rearth * sin(lat_rad)
        θs = range(0, 2π, 100)
        x = [r * cos(θ) for θ in θs]
        y = [r * sin(θ) for θ in θs]
        z_line = fill(z, length(θs))
        p = plot3d!(p, x, y, z_line, lw=0.3, c=:black, alpha=0.3)
    end
	
	return plot!(p, aspect_ratio=:equal, legend=false, colorbar=false, size=(600,600), framestyle=:none, ticks=:none, camera=(0,5))
end 

end