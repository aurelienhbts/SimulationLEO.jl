using LinearAlgebra

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
    T = 2π * sqrt(a^3 / μ)      # Periode orbitale en utilisant la 3e loi de Kepler
    ts = range(0, T; length=n)  # Pour évaluer la couverture à 100 différents moments
    s = 0.0
    for t in ts
        s += coverage_fraction(sats, t, latmin, latmax, eps_deg)
    end
    return s / length(ts)
end

"""
eval_constellation(vec, F, i_deg, a, eps_deg)

Évalue une constellation définie par `vec` et retourne la couverture moyenne sur une période et le nombre total de satellites.
"""
function eval_constellation(vec, F, i_deg, a, eps_deg)
    sats = myconstellation(vec, F, i_deg, a)
    cov = mean_coverage_fraction(sats, -i_deg, i_deg, eps_deg)
    N = length(sats)
    return cov, N
end