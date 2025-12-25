"""
    random_vec(P, N)

Génère un vecteur aléatoire de longueur P représentant une répartition libre
de N satellites parmi P plans orbitaux.

# Arguments
- P : Nombre total de plans orbitaux disponibles.
- N : Nombre total de satellites à répartir.

# Valeur retournée
- Un vecteur de longueur P contenant la distribution aléatoire des N satellites.
"""
function random_vec(P, N)
    v = zeros(Int, P)
    for _ in 1:N
        v[rand(1:P)] += 1
    end
    return v
end

"""
    mutate_vec_fixedN(vec; p_mut=0.3)

Effectue une mutation du vecteur `vec` en déplaçant éventuellement des satellites
d'un plan orbital vers un autre.

# Argument
- vec     : Vecteur de taille P représentant la répartition actuelle des satellites.

# Paramètre optionnel
- p_mut   : Probabilité qu'une mutation soit appliquée à chacun des essais effectués à travers les P plans orbitaux.

# Valeur retournée
- Un nouveau vecteur muté, basé sur `vec`, sans modifier l'original.
"""
function mutate_vec_fixedN(vec; p_mut=0.3)
    P = length(vec)
    v = copy(vec)
    for _ in 1:P
        rand() < p_mut || continue
        i = rand(1:P)
        j = rand(1:P)
        i == j && continue 
        v[i] > 0 || continue # Il faut evidement qu'il y ai un satellite sur le plan i
        v[i] -= 1
        v[j] += 1
    end
    return v
end

const FITCACHE_fixedN = Dict{Tuple{Vararg{Int}}, Float64}()

"""
    fitness_fixedN(vec, F, i_deg, a, eps_deg; nsats=1, Cmin=75.0, Pbonus=true)

Évalue la qualité d'une constellation candidate.

# Arguments
- vec      : Vecteur de taille P indiquant le nombre de satellites dans chaque plan orbital.
- F        : Paramètre de phasage (Walker-Delta) utilisé par `eval_constellation`.
- i_deg    : Inclinaison orbitale en degrés.
- a        : Demi-grand axe de l'orbite (en mètres), généralement Re + altitude.
- eps_deg  : Angle d'élévation minimal (en degrés) nécessaire pour qu'un satellite couvre un point au sol.
- grid_ga  : Stucture `GroundGrid` contenant la grille lat/lon et des valeurs précalculées.

# Paramètres optionnels
- nsats    : Nombre de satellites minimal qui doivent couvrir chaque point au sol.
- Cmin     : Seuil minimal de couverture acceptable. Si la couverture retournée par `eval_constellation` est < Cmin, une pénalité de -100 est appliquée.
- Pbonus   : Si true, récompense légèrement les vecteurs qui utilisent moins de plans orbitaux (bonus proportionnel au nombre de plans vides).

# Valeur retournée
- fit      : Score de qualité de la constellation.  
             Maximisé lorsque la couverture est haute et que la structure utilise peu de plans.
"""
function fitness_fixedN(vec, F, i_deg, a, eps_deg, grid_ga; nsats=1, Cmin=75.0, Pbonus=true)
    
    key = Tuple(vec)
    if haskey(FITCACHE_fixedN, key)
        return FITCACHE_fixedN[key]
    end
    # Evaluation grossière de la configuration (ok pour l'algorithme)
    cov, N = eval_constellation_GA(vec, F, i_deg, a, eps_deg, grid_ga; n=10, nsats=nsats)

    n_used = count(!iszero, vec)
    bonus = Pbonus ? 0.2 * (length(vec) - n_used) : 0 # Si Pbonus=true, on ajoute un bonus en fonction du nombre de plans vides

    fit = cov < Cmin ? cov - 100.0 : cov + bonus # Si cov < Cmin, on applique une forte pénalité
    FITCACHE_fixedN[key] = fit
    return fit
end

"""
    evolve_vec_fixedN(P, N, F, i_deg, a, eps_deg; nsats=1, popsize=20, generations=30, Cmin=0.0, Pbonus=true)

Algorithme génétique simple pour optimiser la répartition de N satellites sur P plans orbitaux.

# Arguments
- P        : Nombre total de plans orbitaux possibles.
- N        : Nombre total de satellites à répartir sur les P plans.
- F        : Paramètre de phasage (Walker-Delta) utilisé dans `eval_constellation` pour définir la géométrie de la constellation.
- i_deg    : Inclinaison orbitale en degrés.
- a        : Demi-grand axe de l'orbite (en mètres), typiquement Re + altitude.
- eps_deg  : Angle d'élévation minimal (en degrés) pour considérer qu'un satellite couvre un point au sol.

# Paramètres optionnels
- nsats        : Nombre de satellites minimal qui doivent couvrir chaque point au sol.
- popsize      : Taille de la population de vecteurs candidats (nombre de solutions évaluées par génération).
- generations  : Nombre de générations de l'algorithme génétique (profondeur de la recherche).
- Cmin      : Seuil minimal de couverture ; si la couverture moyenne est < Cmin, une pénalité forte est appliquée dans `fitness`.
- Pbonus       : Si true, ajoute un bonus dans `fitness` pour les configurations qui utilisent moins de plans orbitaux (favorise les plans vides).

# Valeurs de retour
- best_vec : Vecteur de taille P contenant le nombre de satellites par plan pour la meilleure configuration trouvée.
- cov      : Couverture moyenne correspondante (telle que renvoyée par `eval_constellation`, calculée avec une résolution plus fine à la fin).
"""
function evolve_vec_fixedN(P, N, F, i_deg, a, eps_deg; nsats=1, popsize=20, generations=30, Cmin=75.0, Pbonus=true)

    grid_ga = GroundGrid(-i_deg, i_deg; dlat=6, dlon=6) # Initialisation de GroundGrid

    population = Vector{Vector{Int}}(undef, popsize)
    Threads.@threads for i in 1:popsize
        population[i] = random_vec(P, N) # Population initiale
    end
    best_vec = population[1]
    best_fit = fitness_fixedN(best_vec, F, i_deg, a, eps_deg, grid_ga; nsats=nsats, Cmin=Cmin, Pbonus=Pbonus) # Initialisation de best_fit
    fits = Vector{Float64}(undef, popsize)

    for _ in 1:generations
        Threads.@threads for i in 1:popsize
            fits[i] = fitness_fixedN(population[i], F, i_deg, a, eps_deg, grid_ga; nsats=nsats, Cmin=Cmin, Pbonus=Pbonus)
        end

        order = sortperm(fits, rev=true) # Classement par fitness
        elite = population[order[1:clamp(popsize÷4, 1, popsize)]] # On prend que 1/4 des vecteurs (les meilleurs)
        if fits[order[1]] > best_fit
            best_fit = fits[order[1]] # Mise à jour de best_fit
            best_vec = elite[1] # Mise à jour de best_vec
        end

        newpop = copy(elite)
        while length(newpop) < popsize
            p = elite[rand(1:end)]
            child = mutate_vec_fixedN(p) # Mutation d'un vecteur aléatoire
            push!(newpop, child) # Ajout du mutant
        end
        population = newpop
    end

    cov, _ = eval_constellation(best_vec, F, i_deg, a, eps_deg; n=75, dlat=1, dlon=1, nsats=nsats) # Plus fin pour le cov final
    return best_vec, cov
end
