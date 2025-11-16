using LinearAlgebra

"""
visible(r_ecef, lat_deg, lon_deg, eps_deg)

Vérifie si un satellite est visible depuis un point au sol donné,
en tenant compte d'un angle d'élévation minimal (ε).
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
coverage_fraction(sats, t, latmin, latmax, eps_deg; dlat=2, dlon=2)

Calcule la fraction (%) de points visibles au moins par un satellite **à un instant t donné**.
Tient compte de l'élévation minimale nécessaire pour voir les satellites depuis le sol.
"""
function coverage_fraction(sats, t, latmin, latmax, eps_deg; dlat=2, dlon=2)

    r_ecef = [ecef_from_eci(eci_pos(s, t), t) for s in sats]
    ρs = map(norm, r_ecef)
    cosψmax_s = (Re ./ ρs) .* cos(deg2rad(eps_deg)) # On précalcule cosψmax (angle entre la direction du satellite et celle du point au sol)

    lats = collect(latmin:dlat:latmax)
    lons = collect(-180:dlon:180)

    nthreads = Threads.nthreads()
    covered_threads = fill(0, nthreads)
    pts_threads = fill(0, nthreads)
    
    # Pour des raisons de performances, la fonction 'visible' est directement implémentée ici
    Threads.@threads for i in eachindex(lats)
        tid = Threads.threadid()
        lat = lats[i]
        ϕ = deg2rad(lat)
        cosϕ = cos(ϕ) # On précalcule cosϕ pour ne pas le calculer plusieurs fois
        sinϕ = sin(ϕ) # On précalcule sinϕ pour ne pas le calculer plusieurs fois

        @inbounds for lon in lons
            λ = deg2rad(lon)
            cosλ = cos(λ)
            sinλ = sin(λ)

            # Vecteur qui pointe du centre de la Terre vers le point au sol défini par sa latitude et sa longitude. (repère ECEF)
            gx = cosϕ * cosλ
            gy = cosϕ * sinλ
            gz = sinϕ

            pts_threads[tid] += 1

            for k in eachindex(r_ecef)
                rx, ry, rz = r_ecef[k]
                cosγ = (rx*gx + ry*gy + rz*gz) / ρs[k]
                if cosγ >= cosψmax_s[k]
                    covered_threads[tid] += 1
                    break # Break si il y a dejà 1 satellite au point (pas besoin de calculer plus)
                end
            end
        end
    end

    covered = sum(covered_threads)
    pts = sum(pts_threads)
    return 100 * covered / pts
end

"""
mean_coverage_fraction(sats, latmin, latmax, eps_deg; n=100, dlat=2, dlon=2)

Calcule la fraction (%) des points visibles au moins par un satellite en moyenne sur une période orbitale.
Tient compte de l'élévation minimale nécessaire pour voir les satellites depuis le sol (ε).
"""
function mean_coverage_fraction(sats, latmin, latmax, eps_deg; n=100, dlat=2, dlon=2)
    a = sats[1].a
    T = 2π * sqrt(a^3 / μ) # Période en utilisant la 3e loi de Kepler
    ts = range(0, T; length=n)
    s = 0.0
    for t in ts
        s += coverage_fraction(sats, t, latmin, latmax, eps_deg; dlat=dlat, dlon=dlon)
    end
    return s / length(ts)
end

"""
eval_constellation(vec, F, i_deg, a, eps_deg)

Évalue une constellation définie par `vec` et retourne la couverture moyenne sur une période et le nombre total de satellites.
"""
function eval_constellation(vec, F, i_deg, a, eps_deg; n=100, dlat=2, dlon=2)
    sats = myconstellation(vec, F, i_deg, a)
    cov = mean_coverage_fraction(sats, -i_deg, i_deg, eps_deg; n=n, dlat=dlat, dlon=dlon)
    N = length(sats)
    return cov, N
end

"""
random_vec(P, N)

Génère un vecteur aléatoire de longueur P contenant une répartition libre
de N satellites parmi P plans orbitaux.
"""
function random_vec(P, N)
    vec = zeros(Int, P)              # Compteur de satellites par plan
    for _ in 1:N
        vec[rand(1:P)] += 1          # Ajoute un satellite à un plan choisi au hasard
    end
    return vec
end

"""
mutate_vec(vec; p_mut=0.3)

Effectue des mutations sur `vec` en déplaçant des satellites d'un plan à un autre
selon la probabilité `p_mut`.
"""
function mutate_vec!(vec; p_mut=0.3)
    P = length(vec)
    v = copy(vec)

    for _ in 1:P
        #rand() < p_mut || continue     # Mutation avec probabilité p_mut
        i = rand(1:P)                  # Plan d'où on retire
        j = rand(1:P)                  # Plan où on ajoute
        i == j && continue

        v[i] ≥ 1 || continue # Pour ne pas avoir de plan avec des valeurs négatives
        v[i] -= 1
        v[j] += 1
    end
    return v
end

"""
fitness(vec, F, i_deg, a, eps_deg; Cmin=75.0)

Évalue la qualité d'une constellation.  
- Si la couverture `cov` est sous un seuil `Cmin`, on applique une forte pénalité.  
- Sinon, on maximise `cov` tout en pénalisant faiblement le nombre total de satellites `N`.
"""
function fitness(vec, F, i_deg, a, eps_deg; Cmin=75.0)
    cov, N = eval_constellation(vec, F, i_deg, a, eps_deg; n=20, dlat=4, dlon=4) # Performance actuelle (couverture grossière)

    if cov < Cmin
        return cov - 1000.0    # Forte pénalité si le seuil n'est pas atteint
    else
        return cov # Je pourrais faire 'cov - 0.1*N' mais N es fixe dans l'implémentation actuelle
    end
end

"""
evolve_vec(P, N, F, i_deg, a, eps_deg; popsize=20, generations=30, Cmin=0.0)

Algorithme génétique simple pour optimiser la répartition de N satellites sur P plans orbitaux. Initialise une population aléatoire, sélectionne les meilleures configurations selon `fitness`, les fait muter et retourne le meilleur vecteur trouvé ainsi que sa couverture.
"""
function evolve_vec(P, N, F, i_deg, a, eps_deg; popsize=20, generations=30, Cmin=0.0)
	
    population = [random_vec(P, N) for _ in 1:popsize]              # Population initiale
	
    best_vec = copy(population[1])
    best_fit = fitness(best_vec, F, i_deg, a, eps_deg; Cmin=Cmin)   # Initialisation de best_fit

    for _ in 1:generations
        fits = [fitness(v, F, i_deg, a, eps_deg; Cmin=Cmin) for v in population]
        order = sortperm(fits, rev=true)                            # Classement par fitness
        population = population[order]
        fits = fits[order]

        if fits[1] > best_fit
            best_fit = fits[1]
            best_vec = copy(population[1])                          # Mise à jour de best_fit
        end

        elite = population[1:clamp(popsize ÷ 4, 1, popsize)]        # Élites conservées
        newpop = copy(elite)
        while length(newpop) < popsize
            parent = elite[rand(1:length(elite))]
            child = mutate_vec!(parent)                              # Mutation d'un parent
            push!(newpop, child)
        end
        population = newpop
    end

    cov, _ = eval_constellation(best_vec, F, i_deg, a, eps_deg)
    return best_vec, cov
end