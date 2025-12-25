"""
  GroundGrid

Structure contenant une grille lat/lon ainsi que les valeurs trigonométriques associées aux latitudes et longitudes.
Elle est utilisée pour accélérer les calculs de couverture en évitant le recalcul répété des fonctions trigonométriques lors des évaluations successives de constellations.

# Champs
- lats  : Vecteur des latitudes échantillonnées de la grille sol (en degrés).
- lons  : Vecteur des longitudes échantillonnées de la grille sol (en degrés).
- cosϕ  : Valeurs de cosinus des latitudes correspondantes.
- sinϕ  : Valeurs de sinus des latitudes correspondantes.
- cosλ  : Valeurs de cosinus des longitudes correspondantes.
- sinλ  : Valeurs de sinus des longitudes correspondantes.
"""
struct GroundGrid
    lats::Vector{Float64}
    lons::Vector{Float64}
    cosϕ::Vector{Float64}
    sinϕ::Vector{Float64}
    cosλ::Vector{Float64}
    sinλ::Vector{Float64}
end

"""
    GroundGrid(latmin, latmax; dlat=2, dlon=2)

Construit une grille sol régulière comprise entre `latmin` et `latmax`, et précalcule les valeurs trigonométriques associées aux latitudes et longitudes de la grille.

# Arguments
- latmin   : Latitude minimale des points au sol à tester (en degrés).
- latmax   : Latitude maximale des points au sol à tester (en degrés).

# Paramètres optionnels
- dlat     : Pas d'échantillonnage en latitude (en degrés).
- dlon     : Pas d'échantillonnage en longitude (en degrés).

# Valeur retournée
- Une structure `GroundGrid` contenant la grille sol et les valeurs trigonométriques précalculées associées.
"""
function GroundGrid(latmin, latmax; dlat=2, dlon=2)
    lats = collect(float(latmin):float(dlat):float(latmax))
    lons = collect(-180.0:float(dlon):180.0)

    cosϕ = Vector{Float64}(undef, length(lats))
    sinϕ = Vector{Float64}(undef, length(lats))
    for i in eachindex(lats)
        ϕ = deg2rad(lats[i])
        cosϕ[i] = cos(ϕ)
        sinϕ[i] = sin(ϕ)
    end

    cosλ = Vector{Float64}(undef, length(lons))
    sinλ = Vector{Float64}(undef, length(lons))
    for j in eachindex(lons)
        λ = deg2rad(lons[j])
        cosλ[j] = cos(λ)
        sinλ[j] = sin(λ)
    end

    return GroundGrid(lats, lons, cosϕ, sinϕ, cosλ, sinλ)
end

"""
    coverage_fraction_GA(sats, t, grid::GroundGrid, eps_deg; nsats=1)

Version modifiée de `coverage_fraction`.
Calcule la fraction (%) de points visibles par nsats satellites **à l'instant t**.
La visibilité tient compte de l'élévation minimale `eps_deg`, c'est-à-dire de l'angle sous lequel un point au sol doit voir le satellite pour être considéré comme couvert.

# Arguments
- sats     : Liste des satellites constituant la constellation.
- t        : Instant d'évaluation (en secondes) dans le repère inertiel.
- grid_ga  : Stucture `GroundGrid` contenant la grille lat/lon et des valeurs précalculées.
- eps_deg  : Angle d'élévation minimal (en degrés) pour considérer un satellite visible.

# Paramètres optionnels
- nsats    : Nombre de satellites minimal qui doivent couvrir chaque point au sol.

# Valeur retournée
- Pourcentage de points visibles au moins par un satellite à l'instant t.
"""
function coverage_fraction_GA(sats, t, grid::GroundGrid, eps_deg; nsats=1)
    r_ecef = [ecef_from_eci(eci_pos(s, t), t) for s in sats]
    ρs = map(norm, r_ecef)
    cosψmax_s = (Re ./ ρs) .* cos(deg2rad(eps_deg))

    nlon = length(grid.lons)
    covered_thread = zeros(Int, Threads.nthreads())

    Threads.@threads for i in eachindex(grid.lats)
        tid = Threads.threadid()
        cϕ = grid.cosϕ[i]
        sϕ = grid.sinϕ[i]
        local_cov = 0

        @inbounds for j in 1:nlon
            gx = cϕ * grid.cosλ[j]
            gy = cϕ * grid.sinλ[j]
            gz = sϕ

            hit = 0
            @inbounds for k in eachindex(r_ecef)
                rx, ry, rz = r_ecef[k]
                cosγ = (rx*gx + ry*gy + rz*gz) / ρs[k]
                if cosγ >= cosψmax_s[k]
                    hit += 1
                    if hit == nsats
                        local_cov += 1
                        break
                    end
                end
            end
        end

        covered_thread[tid] += local_cov
    end

    covered = sum(covered_thread)
    pts = length(grid.lats) * length(grid.lons)
    return 100 * covered / pts
end

"""
    function mean_coverage_fraction_GA(sats, grid::GroundGrid, eps_deg; n=100, nsats=1)

Version modifiée de `mean_coverage_fraction`.
Calcule la fraction (%) des points visibles par au moins nsats satellites en moyenne sur une période orbitale complète.  
Tient compte de l'élévation minimale `eps_deg` requise pour considérer qu'un point au sol est effectivement couvert.

# Arguments
- sats     : Liste d'objets satellites constituant la constellation.
- grid_ga  : Stucture `GroundGrid` contenant la grille lat/lon et des valeurs précalculées.
- eps_deg  : Angle d'élévation minimal pour considérer qu'un satellite couvre un point.

# Paramètres optionnels
- n        : Nombre de pas de temps pour évaluer une période orbitale.
- nsats    : Nombre de satellites minimal qui doivent couvrir chaque point au sol.

# Valeur retournée
- La couverture moyenne, exprimée en fraction, obtenue en évaluant la couverture
  à n instants répartis sur une période orbitale.
"""
function mean_coverage_fraction_GA(sats, grid::GroundGrid, eps_deg; n=100, nsats=1)
    a = sats[1].a
    T = 2π * sqrt(a^3 / μ)
    ts = range(0, T; length=n)
    s = 0.0
    for t in ts
        s += coverage_fraction_GA(sats, t, grid, eps_deg; nsats=nsats)
    end
    return s / length(ts)
end

"""
    eval_constellation_GA(vec, F, i_deg, a, eps_deg, grid::GroundGrid; n=100, nsats=1)

Version modifiée de `eval_constellation`.
Évalue une constellation décrite par le vecteur `vec` et retournela couverture moyenne obtenue sur une période ainsi que le nombre total de satellites.

# Arguments
- vec     : Vecteur indiquant le nombre de satellites dans chaque plan orbital.
- F       : Paramètre de phasage (Walker-Delta) utilisé pour construire la constellation.
- i_deg   : Inclinaison orbitale en degrés.
- a       : Demi-grand axe de l'orbite (en mètres), généralement Re + altitude.
- eps_deg : Angle d'élévation minimal pour considérer qu'un point au sol est couvert.
- grid_ga  : Stucture `GroundGrid` contenant la grille lat/lon et des valeurs précalculées.

# Paramètres optionnels
- n        : Nombre de pas de temps pour évaluer une période orbitale.
- nsats    : Nombre de satellites minimal qui doivent couvrir chaque point au sol.

# Valeurs retournées
- cov     : Couverture moyenne en pourcentage (fraction de points visibles en moyenne).
- N       : Nombre total de satellites générés à partir de `vec`.
"""
function eval_constellation_GA(vec, F, i_deg, a, eps_deg, grid::GroundGrid; n=100, nsats=1)
    sats = myconstellation(vec, F, i_deg, a)
    cov = mean_coverage_fraction_GA(sats, grid, eps_deg; n=n, nsats=nsats)
    return cov, length(sats)
end