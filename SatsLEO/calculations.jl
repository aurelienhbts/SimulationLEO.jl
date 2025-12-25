"""
    visible(r_ecef, lat_deg, lon_deg, eps_deg)

Détermine si un satellite est visible depuis un point donné au sol, en se basant sur l'angle d'élévation minimal `eps_deg`.  
La fonction calcule l'angle entre le vecteur position du satellite (en ECEF) et le vecteur normal au point au sol.  
Si cet angle est compatible avec l'élévation minimale requise, le satellite est considéré comme visible.

# Arguments
- r_ecef   : Vecteur position du satellite dans le repère ECEF.
- lat_deg  : Latitude du point au sol (en degrés).
- lon_deg  : Longitude du point au sol (en degrés).
- eps_deg  : Angle d'élévation minimal (en degrés) pour considérer le satellite visible.

# Valeur retournée
- true si le satellite satisfait la contrainte d'élévation minimale.
- false sinon.
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
    coverage_fraction(sats, t, latmin, latmax, eps_deg; dlat=2, dlon=2, nsats=1)

Calcule la fraction (%) de points visibles par nsats satellites **à l'instant t**.
La visibilité tient compte de l'élévation minimale `eps_deg`, c'est-à-dire de l'angle sous lequel un point au sol doit voir le satellite pour être considéré comme couvert.

# Arguments
- sats     : Liste des satellites constituant la constellation.
- t        : Instant d'évaluation (en secondes) dans le repère inertiel.
- latmin   : Latitude minimale des points au sol à tester.
- latmax   : Latitude maximale des points au sol à tester.
- eps_deg  : Angle d'élévation minimal (en degrés) pour considérer un satellite visible.
- nsats    : Nombre de satellites minimal qui doivent couvrir chaque point au sol.

# Paramètres optionnels
- dlat     : Pas d'échantillonnage en latitude (en degrés).
- dlon     : Pas d'échantillonnage en longitude (en degrés).

# Valeur retournée
- Pourcentage de points visibles au moins par un satellite à l'instant t.
"""
function coverage_fraction(sats, t, latmin, latmax, eps_deg; dlat=2, dlon=2, nsats=1)
    r_ecef = [ecef_from_eci(eci_pos(s, t), t) for s in sats]
    ρs = map(norm, r_ecef)
    cosψmax_s = (Re ./ ρs) .* cos(deg2rad(eps_deg)) # On précalcule cosψmax (angle entre la direction du satellite et celle du point au sol)

    lats = collect(latmin:dlat:latmax)
    lons = collect(-180:dlon:180)

    cosϕ = similar(lats, Float64)
    sinϕ = similar(lats, Float64)
    for i in eachindex(lats)
        ϕ = deg2rad(lats[i])
        cosϕ[i] = cos(ϕ) # On précalcule cosϕ pour ne pas le calculer plusieurs fois
        sinϕ[i] = sin(ϕ) # On précalcule sinϕ pour ne pas le calculer plusieurs fois
    end

    cosλ = similar(lons, Float64)
    sinλ = similar(lons, Float64)
    for j in eachindex(lons)
        λ = deg2rad(lons[j])
        cosλ[j] = cos(λ) # On précalcule cosλ pour ne pas le calculer plusieurs fois
        sinλ[j] = sin(λ) # On précalcule sinλ pour ne pas le calculer plusieurs fois
    end

    nlon = length(lons)
    covered_thread = zeros(Int, Threads.nthreads())

    Threads.@threads for i in eachindex(lats)
        tid = Threads.threadid()
        cϕ = cosϕ[i]
        sϕ = sinϕ[i]
        local_cov = 0

        @inbounds for j in 1:nlon
            # Vecteur qui pointe du centre de la Terre vers le point au sol défini par sa latitude et sa longitude. (repère ECEF)
            gx = cϕ * cosλ[j]
            gy = cϕ * sinλ[j]
            gz = sϕ

            hit = 0
            @inbounds for k in eachindex(r_ecef)
                rx, ry, rz = r_ecef[k]
                cosγ = (rx*gx + ry*gy + rz*gz) / ρs[k]
                if cosγ >= cosψmax_s[k]
                    hit += 1
                    if hit == nsats # On compte le point si il y a plus de nsat satellites qui le couvrent
                        local_cov += 1
                        break
                    end
                end
            end
        end

        covered_thread[tid] += local_cov
    end

    covered = sum(covered_thread)
    pts = length(lats) * length(lons)
    return 100 * covered / pts
end

"""
    mean_coverage_fraction(sats, latmin, latmax, eps_deg; n=100, dlat=2, dlon=2, nsats=1)

Calcule la fraction (%) des points visibles par au moins nsats satellites en moyenne sur une période orbitale complète.  
Tient compte de l'élévation minimale `eps_deg` requise pour considérer qu'un point au sol est effectivement couvert.

# Arguments
- sats     : Liste d'objets satellites constituant la constellation.
- latmin   : Latitude minimale (en degrés) des points au sol à évaluer.
- latmax   : Latitude maximale (en degrés) des points au sol à évaluer.
- eps_deg  : Angle d'élévation minimal pour considérer qu'un satellite couvre un point.

# Paramètres optionnels
- n        : Nombre de pas de temps pour évaluer une période orbitale.
- dlat     : Résolution en latitude (en degrés) pour l'échantillonnage au sol.
- dlon     : Résolution en longitude (en degrés) pour l'échantillonnage au sol.
- nsats    : Nombre de satellites minimal qui doivent couvrir chaque point au sol.

# Valeur retournée
- La couverture moyenne, exprimée en fraction, obtenue en évaluant la couverture
  à n instants répartis sur une période orbitale.
"""
function mean_coverage_fraction(sats, latmin, latmax, eps_deg; n=100, dlat=2, dlon=2, nsats=1)
    a = sats[1].a
    T = 2π * sqrt(a^3 / μ) # Période en utilisant la 3e loi de Kepler
    ts = range(0, T; length=n)
    s = 0.0
    for t in ts
        s += coverage_fraction(sats, t, latmin, latmax, eps_deg; dlat=dlat, dlon=dlon, nsats=nsats)
    end
    return s / length(ts)
end

"""
    eval_constellation(vec, F, i_deg, a, eps_deg; n=100, dlat=2, dlon=2, nsats=1)

Évalue une constellation décrite par le vecteur `vec` et retourne la couverture moyenne obtenue sur une période ainsi que le nombre total de satellites.

# Arguments
- vec     : Vecteur indiquant le nombre de satellites dans chaque plan orbital.
- F       : Paramètre de phasage (Walker-Delta) utilisé pour construire la constellation.
- i_deg   : Inclinaison orbitale en degrés.
- a       : Demi-grand axe de l'orbite (en mètres), généralement Re + altitude.
- eps_deg : Angle d'élévation minimal pour considérer qu'un point au sol est couvert.

# Paramètres optionnels
- n        : Nombre de pas de temps pour évaluer une période orbitale.
- dlat    : Résolution en latitude (en degrés) pour les points tests au sol.
- dlon    : Résolution en longitude (en degrés) pour les points tests au sol.
- nsats    : Nombre de satellites minimal qui doivent couvrir chaque point au sol.

# Valeurs retournées
- cov     : Couverture moyenne en pourcentage (fraction de points visibles en moyenne).
- N       : Nombre total de satellites générés à partir de `vec`.
"""
function eval_constellation(vec, F, i_deg, a, eps_deg; n=100, dlat=2, dlon=2, nsats=1)
    sats = myconstellation(vec, F, i_deg, a)
    cov = mean_coverage_fraction(sats, -i_deg, i_deg, eps_deg; n=n, dlat=dlat, dlon=dlon, nsats=nsats)
    N = length(sats)
    return cov, N
end