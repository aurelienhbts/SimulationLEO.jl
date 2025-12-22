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

const FITCACHE_GLOBAL = Dict{Tuple{Vararg{Int}}, Float64}()
const FITLOCK = Threads.SpinLock()
const FITCACHE_LOCAL = [Dict{Tuple{Vararg{Int}}, Float64}() for _ in 1:Threads.nthreads()]

function flush_fitcache!()
    Threads.lock(FITLOCK)
    for d in FITCACHE_LOCAL
        for (k, v) in d
            FITCACHE_GLOBAL[k] = v
        end
        empty!(d)
    end
    Threads.unlock(FITLOCK)
    return nothing
end

"""
    fitness_fixedN(vec, F, i_deg, a, eps_deg; nsats=1, Cmin=75.0, Pbonus=true)

Évalue la qualité d'une constellation candidate.

# Arguments
- vec      : Vecteur de taille P indiquant le nombre de satellites dans chaque plan orbital.
- F        : Paramètre de phasage (Walker-Delta) utilisé par `eval_constellation`.
- i_deg    : Inclinaison orbitale en degrés.
- a        : Demi-grand axe de l'orbite (en mètres), généralement Re + altitude.
- eps_deg  : Angle d'élévation minimal (en degrés) nécessaire pour qu'un satellite couvre un point au sol.

# Paramètres optionnels
- nsats    : Nombre de satellites minimal qui doivent couvrir chaque point au sol.
- Cmin     : Seuil minimal de couverture acceptable. Si la couverture retournée par `eval_constellation` est < Cmin, une pénalité de -100 est appliquée.
- Pbonus   : Si true, récompense légèrement les vecteurs qui utilisent moins de plans orbitaux (bonus proportionnel au nombre de plans vides).

# Valeur retournée
- fit      : Score de qualité de la constellation.  
             Maximisé lorsque la couverture est haute et que la structure utilise peu de plans.
"""
function fitness_fixedN(vec, F, i_deg, a, eps_deg; nsats=1, Cmin=75.0, Pbonus=true)
    key = Tuple(vec)
    tid = Threads.threadid()

    vloc = get(FITCACHE_LOCAL[tid], key, NaN)
    if vloc === vloc
        return vloc
    end

    Threads.lock(FITLOCK)
    vglob = get(FITCACHE_GLOBAL, key, NaN)
    Threads.unlock(FITLOCK)
    if vglob === vglob
        FITCACHE_LOCAL[tid][key] = vglob
        return vglob
    end

    cov, N = eval_constellation(vec, F, i_deg, a, eps_deg; n=10, dlat=6, dlon=6, nsats=nsats)

    n_used = count(!iszero, vec)
    bonus = Pbonus ? 0.2 * (length(vec) - n_used) : 0.0

    fit = cov < Cmin ? cov - 100.0 : cov + bonus
    FITCACHE_LOCAL[tid][key] = fit
    return fit
end

function fitness_fixedN_v1(vec, F, i_deg, a, eps_deg; nsats=1, Cmin=75.0, Pbonus=true)
    
    key = Tuple(vec)
    if haskey(FITCACHE, key)
        return FITCACHE[key]
    end
    # Evaluation grossière de la configuration (ok pour l'algorithme)
    cov, N = eval_constellation(vec, F, i_deg, a, eps_deg; n=10, dlat=6, dlon=6, nsats=nsats)

    n_used = count(!iszero, vec)
    bonus = Pbonus ? 0.2 * (length(vec) - n_used) : 0 # Si Pbonus=true, on ajoute un bonus en fonction du nombre de plans vides

    fit = cov < Cmin ? cov - 100.0 : cov + bonus # Si cov < Cmin, on applique une forte pénalité
    FITCACHE[key] = fit
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

    population = Vector{Vector{Int}}(undef, popsize)
    Threads.@threads for i in 1:popsize
        population[i] = random_vec(P, N) # Population initiale
    end
    best_vec = population[1]
    best_fit = fitness_fixedN(best_vec, F, i_deg, a, eps_deg; nsats=nsats, Cmin=Cmin, Pbonus=Pbonus) # Initialisation de best_fit
    fits = Vector{Float64}(undef, popsize)

    for _ in 1:generations
        Threads.@threads for i in 1:popsize
            fits[i] = fitness_fixedN(population[i], F, i_deg, a, eps_deg; nsats=nsats, Cmin=Cmin, Pbonus=Pbonus)
        end
        flush_fitcache!() # Recombinaison des caches par thread

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


"""
    mutate_vec(vec; p_move=0.4, p_add=0.1, p_rem=0.05)

Applique une mutation à une configuration orbitale.
La mutation peut déplacer, ajouter ou retirer des satellites, permettant ainsi d'explorer des constellations de tailles et de répartitions variées.

# Arguments
- vec      : Vecteur de taille P indiquant le nombre de satellites dans chaque plan orbital.

# Paramètres optionnels
- p_move   : Probabilité de déplacer un satellite d'un plan vers un autre (répété P fois).
- p_add    : Probabilité d'ajouter un satellite dans un plan choisi aléatoirement.
- p_rem    : Probabilité d'enlever un satellite dans un plan non vide.

# Valeur retournée
- v        : Nouveau vecteur muté, basé sur `vec`.
"""
function mutate_vec(vec; p_move=0.4, p_add=0.1, p_rem=0.05)
    P = length(vec)
    v = copy(vec)

    for _ in 1:P
        if rand() < p_move && sum(v) > 1 
            inds = findall(>(0), v) # Plans avec des satellites
            i = rand(inds)
            j = rand(1:P)
            if i != j
                v[i] -= 1
                v[j] += 1
            end
        end
    end

    if rand() < p_add 
        k = rand(1:P)
        v[k] += 1
    end

    if rand() < p_rem && sum(v) > 1 
        inds = findall(>(0), v)
        k = rand(inds)
        v[k] -= 1
    end

    return v
end

"""
    fitness(vec, F, i_deg, a, eps_deg; nsats=1, 
            Ncoef=0.75, Pcoef=0.3, Ctarget=95.0, K=5.0)

Évalue la qualité d'une configuration orbitale variable en nombre de satellites.
Maximisé lorsque la couverture est élevée, que le nombre de satellites reste faible, que l'usage des plans est concentré, et que les constellations “presque saturées” (cov ≈ 100 %) ne gonflent pas artificiellement N.

# Arguments
- vec      : Vecteur de taille P indiquant le nombre de satellites dans chaque plan orbital.
- F        : Paramètre de phasage (Walker-Delta) utilisé par `myconstellation`.
- i_deg    : Inclinaison orbitale en degrés.
- a        : Demi-grand axe de l'orbite (mètres).
- eps_deg  : Angle d'élévation minimal requis pour la visibilité.

# Paramètres optionnels
- nsats    : Nombre de satellites minimal qui doivent couvrir chaque point au sol.
- Ncoef    : Coefficient contrôlant la pénalité liée au nombre total de satellites. Le malus est faible pour N < 17, élevé pour N > 23, et interpolé linéairement entre les deux.
- Pcoef    : Bonus attribué aux configurations comportant des plans orbitaux vides.
- Ctarget  : Couverture minimale souhaitée. Si la couverture retournée est < Ctarget, une pénalité proportionnelle au déficit de couverture est appliquée, sans empêcher l'exploration de bonnes solutions légèrement en dessous.
- K        : Intensité de la pénalité lorsque cov < Ctarget.

# Valeur retournée
- fit      : Score de qualité de la constellation.
"""
function fitness(vec, F, i_deg, a, eps_deg; nsats=1, Ncoef=0.75, Pcoef=0.3, Ctarget=95.0, K=5.0)
    key = Tuple(vec)
    tid = Threads.threadid()

    vloc = get(FITCACHE_LOCAL[tid], key, NaN)
    if vloc === vloc
        return vloc
    end

    Threads.lock(FITLOCK)
    vglob = get(FITCACHE_GLOBAL, key, NaN)
    Threads.unlock(FITLOCK)
    if vglob === vglob
        FITCACHE_LOCAL[tid][key] = vglob
        return vglob
    end

    cov, N = eval_constellation(vec, F, i_deg, a, eps_deg; n=10, dlat=6, dlon=6, nsats=nsats)
    P_empty = length(vec) - count(!iszero, vec)

    if N < 17
        Nmalus = 0.3 * Ncoef
    elseif N > 23
        Nmalus = Ncoef
    else
        t = (N - 17) / 6
        Nmalus = (0.3 + 0.7 * t) * Ncoef
    end

    if cov >= 99.9
        Nmalus *= 3.0
    elseif cov >= 99.5
        Nmalus *= 2.0
    end

    penalty_cov = max(0.0, Ctarget - cov) * K
    fit = cov - Nmalus * N + Pcoef * P_empty - penalty_cov

    FITCACHE_LOCAL[tid][key] = fit
    return fit
end

function fitness_v1(vec, F, i_deg, a, eps_deg; nsats=1, Ncoef=0.75, Pcoef=0.3, Ctarget=95.0, K=5.0)

    key = Tuple(vec)
    if haskey(FITCACHE, key)
        return FITCACHE[key]
    end
    # Evaluation grossière de la configuration (ok pour l'algorithme)
    cov, N = eval_constellation(vec, F, i_deg, a, eps_deg; n=10, dlat=6, dlon=6, nsats=nsats)
    P_empty = length(vec) - count(!iszero, vec) # Nombre de plans vides

    if N < 17 # Faible malus si N est petit
        Nmalus = 0.3 * Ncoef
    elseif N > 23 # Grand malus si N est grand
        Nmalus = Ncoef
    else # interpolation entre 17 et 23
        t = (N - 17) / 6
        Nmalus = (0.3 + 0.7 * t) * Ncoef
    end

    if cov >= 99.9
        Nmalus *= 3.0
    elseif cov >= 99.5
        Nmalus *= 2.0
    end
    penalty_cov = max(0.0, Ctarget - cov) * K

    fit = cov - Nmalus * N + Pcoef * P_empty - penalty_cov
    FITCACHE[key] = fit
    return fit
end

"""
    evolve_vec(P, N_init, F, i_deg, a, eps_deg; nsats=1, popsize=30, generations=40,
               Ncoef=0.75, Pcoef=0.3, Ctarget=95.0, K=5.0, p_move=0.4, p_add=0.1, p_rem=0.05)

Algorithme génétique qui optimise une constellation LEO en laissant varier à la fois
la répartition des satellites par plan et le nombre total de satellites.

# Arguments
- P          : Nombre total de plans orbitaux disponibles.
- N_init     : Nombre total de satellites utilisés pour initialiser la population.
- F          : Paramètre de phasage (Walker-Delta) utilisé par `eval_constellation`.
- i_deg      : Inclinaison orbitale en degrés.
- a          : Demi-grand axe de l'orbite (en mètres), typiquement Re + altitude.
- eps_deg    : Angle d'élévation minimal requis pour la visibilité.

# Paramètres optionnels
- nsats    : Nombre de satellites minimal qui doivent couvrir chaque point au sol.
- popsize    : Taille de la population évoluée à chaque génération.
- generations: Nombre de générations de l'algorithme génétique.
- Ncoef      : Coefficient contrôlant la pénalité liée au nombre total de satellites dans `fitness`.
- Pcoef      : Coefficient contrôlant le bonus associé aux plans orbitaux vides dans `fitness`.
- Ctarget    : Couverture minimale souhaitée utilisée par `fitness` pour pénaliser les constellations trop faibles.
- K          : Intensité de la pénalité lorsque cov < Ctarget.
- p_move     : Probabilité de déplacer un satellite d'un plan vers un autre lors de `mutate_vec`.
- p_add      : Probabilité d'ajouter un satellite dans un plan lors de `mutate_vec`.
- p_rem      : Probabilité d'enlever un satellite dans un plan non vide lors de `mutate_vec`.

# Valeurs retournées
- best_vec   : Vecteur (taille P) décrivant la meilleure répartition trouvée.
- cov_final  : Couverture moyenne finale associée à `best_vec` (évaluée finement).
- N_final    : Nombre total de satellites de la configuration optimale retenue.
"""
function evolve_vec(P, N_init, F, i_deg, a, eps_deg; nsats=1, popsize=30, generations=40, 
                    Ncoef=0.75, Pcoef=0.3, Ctarget=95.0, K=5.0, p_move=0.4, p_add=0.1, p_rem=0.05)

    population = Vector{Vector{Int}}(undef, popsize)
    Threads.@threads for i in 1:popsize
        population[i] = random_vec(P, N_init) # Population initiale
    end
    best_vec = population[1] # Initialisation de best_vec
    best_fit = fitness(best_vec, F, i_deg, a, eps_deg; nsats=nsats, Ncoef=Ncoef, Pcoef=Pcoef, Ctarget=Ctarget, K=K) # Initialisation de best_fit (pour best_vec)
    fits = Vector{Float64}(undef, popsize)

    for _ in 1:generations
        Threads.@threads for i in 1:popsize
            fits[i] = fitness(population[i], F, i_deg, a, eps_deg; nsats=nsats, Ncoef=Ncoef, Pcoef=Pcoef, Ctarget=Ctarget)
        end
        flush_fitcache!() # Recombinaison des caches par thread

        order = sortperm(fits, rev=true) # On ordonne la population en fonction des fitness
        elite_count = clamp(popsize ÷ 4, 1, popsize)
        elite = population[order[1:elite_count]] # On prend que les meilleurs

        if fits[order[1]] > best_fit
            best_fit = fits[order[1]]
            best_vec = elite[1]
        end

        newpop = Vector{typeof(best_vec)}()
        append!(newpop, elite) # Ajout des elites dans la nouvelle population

        while length(newpop) < popsize
            parent = elite[rand(1:elite_count)] # Mutation d'un parent aléatoire
            child = mutate_vec(parent; p_move=p_move, p_add=p_add, p_rem=p_rem)
            push!(newpop, child) # Ajout de la nouvelle configuration
        end

        for _ in 1:clamp(popsize ÷ 10, 1, popsize - length(newpop))
            push!(newpop, random_vec(P, N_init)) # Ajout de vecteurs aléatoires 
        end

        population = newpop
    end

    cov_final, N_final = eval_constellation(best_vec, F, i_deg, a, eps_deg; n=75, dlat=1, dlon=1, nsats=nsats)
    return best_vec, cov_final, N_final
end