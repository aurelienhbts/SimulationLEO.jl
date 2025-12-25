"""
    mutate_vec(vec; p_move=0.4, p_add=0.1, p_rem=0.05)

Applique une mutation à une configuration orbitale.
La mutation peut déplacer, ajouter ou retirer des satellites, permettant ainsi d'explorer des constellations de tailles et de répartitions variées.

# Arguments
- vec      : Vecteur de taille P indiquant le nombre de satellites dans chaque plan orbital.

# Paramètres optionnels
- p_move     : Probabilité de déplacer un satellite d'un plan vers un autre (répété P fois).
- p_add_max  : Probabilité d'ajouter un satellite dans un plan choisi aléatoirement.
- p_rem_max  : Probabilité d'enlever un satellite dans un plan non vide.
- Nmax       : Nombre maximal de satellites souhaités.

# Valeur retournée
- v        : Nouveau vecteur muté, basé sur `vec`.
"""
function mutate_vec(vec, best_cov, Ctarget; p_move=0.4, p_add_max=0.3, p_rem_max=0.2, Nmax=30)
    P = length(vec)
    v = copy(vec)
    N = sum(v)

    gap = (best_cov - Ctarget) / Ctarget # p_add et p_rem variables en fonction du gap
    p_add = p_add_max * clamp(-gap, 0.0, 1.0)
    p_rem = p_rem_max * clamp(gap, 0.0, 1.0)

    for _ in 1:P
        if rand() < p_move && N > 1
            i = rand(1:P)
            tries = 0
            while v[i] == 0 && tries < P # On prend que les plans vides
                i = rand(1:P)
                tries += 1
            end
            v[i] == 0 && continue # On continue si on a pas trouvé de plan non-vide
            j = rand(1:P)
            i == j && continue
            v[i] -= 1 # Retrait
            v[j] += 1 # Ajout
        end
    end

    if rand() < p_add && N < Nmax # Ajout d'un satellite sur un plan au hasard
        k = rand(1:P)
        v[k] += 1
        N += 1
    end

    if rand() < p_rem && N > 1 # Retrait d'un satellite sur un plan au hasard
        k = rand(1:P)
        tries = 0
        while v[k] == 0 && tries < P
            k = rand(1:P)
            tries += 1
        end
        v[k] == 0 && return v # On ne retire pas si on a pas trouvé de plan non-vide
        v[k] -= 1
    end

    return v
end

const FITCACHE = Dict{Tuple{Vararg{Int}}, Tuple{Float64,Float64}}()

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
- grid_ga  : Stucture `GroundGrid` contenant la grille lat/lon et des valeurs précalculées.

# Paramètres optionnels
- nsats    : Nombre de satellites minimal qui doivent couvrir chaque point au sol.
- Ncoef    : Coefficient contrôlant la pénalité liée au nombre total de satellites. Le malus est faible pour N < 17, élevé pour N > 23, et interpolé linéairement entre les deux.
- Pcoef    : Bonus attribué aux configurations comportant des plans orbitaux vides.
- Ctarget  : Couverture minimale souhaitée. Si la couverture retournée est < Ctarget, une pénalité proportionnelle au déficit de couverture est appliquée, sans empêcher l'exploration de bonnes solutions légèrement en dessous.
- K        : Intensité de la pénalité lorsque cov < Ctarget.

# Valeur retournée
- cov      : Coverage de la constellation.
- fit      : Score de qualité de la constellation.
"""
function fitness(vec, F, i_deg, a, eps_deg, grid_ga; nsats=1, Ncoef=0.75, Pcoef=0.3, Ctarget=95.0, K=5.0)

    key = Tuple(vec)
    if haskey(FITCACHE, key)
        return FITCACHE[key]
    end
    # Evaluation grossière de la configuration (ok pour l'algorithme)
    cov, N = eval_constellation_GA(vec, F, i_deg, a, eps_deg, grid_ga; n=10, nsats=nsats)
    P_empty = length(vec) - count(!iszero, vec) # Nombre de plans vides

    if N < 15 # Faible malus si N est petit
        Nmalus = 0.3 * Ncoef
    elseif N > 30 # Grand malus si N est grand
        Nmalus = Ncoef
    else # interpolation entre 17 et 23
        t = (N - 15) / 15
        Nmalus = (0.3 + 0.7 * t) * Ncoef
    end

    if cov >= 99.9999 # Augmenter le malus si la cov depasse 100%
        Nmalus *= 3.0
    end
    penalty_cov = max(0.0, Ctarget - cov) * K

    fit = cov - Nmalus * N + Pcoef * P_empty - penalty_cov
    FITCACHE[key] = cov, fit
    return cov, fit
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
- Nmax       : Nombre maximal de satellites souhaités.
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
function evolve_vec(P, N_init, F, i_deg, a, eps_deg; nsats=1, popsize=30, generations=40, Nmax=30,
                    Ncoef=0.75, Pcoef=0.3, Ctarget=95.0, K=5.0, p_move=0.4, p_add_max=0.3, p_rem_max=0.2)

    grid_ga = GroundGrid(-i_deg, i_deg; dlat=6, dlon=6) # Initialisation de GroundGrid

    population = Vector{Vector{Int}}(undef, popsize)
    Threads.@threads for i in 1:popsize
        population[i] = random_vec(P, N_init) # Population initiale
    end

    best_vec = population[1] # Initialisation de best_vec
    best_cov, best_fit = fitness(best_vec, F, i_deg, a, eps_deg, grid_ga; nsats=nsats, Ncoef=Ncoef, Pcoef=Pcoef, Ctarget=Ctarget, K=K)

    fits = Vector{Float64}(undef, popsize) # Vecteur avec les valeurs de fit de la population actuelle
    covs = Vector{Float64}(undef, popsize) # Vecteur avec les valeurs de cov de la population actuelle

    for _ in 1:generations
        Threads.@threads for i in 1:popsize
            cov, fit = fitness(population[i], F, i_deg, a, eps_deg, grid_ga; nsats=nsats, Ncoef=Ncoef, Pcoef=Pcoef, Ctarget=Ctarget, K=K)
            covs[i] = cov
            fits[i] = fit
        end

        order = sortperm(fits, rev=true)  # On ordonne la population en fonction des fitness
        elite_count = clamp(popsize ÷ 4, 1, popsize)
        elite = population[order[1:elite_count]] # On prend que les meilleurs

        if fits[order[1]] > best_fit
            best_fit = fits[order[1]]
            best_cov = covs[order[1]]
            best_vec = elite[1]
        end

        newpop = Vector{typeof(best_vec)}()
        append!(newpop, elite) # Ajout des elites dans la nouvelle population

        while length(newpop) < popsize - 1
            parent = elite[rand(1:elite_count)] # Mutation d'un parent aléatoire
            child = mutate_vec(parent, best_cov, Ctarget; p_move=p_move, p_add_max=p_add_max, p_rem_max=p_rem_max, Nmax=Nmax)
            push!(newpop, child) # Ajout de la nouvelle configuration
        end

        push!(newpop, random_vec(P, N_init)) # Ajout d'un vecteur aléatoire (avec un peu de chance, on découvre une nouvelle partie de l'espace des configurations)
        population = newpop
    end

    cov_final, N_final = eval_constellation(best_vec, F, i_deg, a, eps_deg; n=75, dlat=1, dlon=1, nsats=nsats)
    return best_vec, cov_final, N_final
end