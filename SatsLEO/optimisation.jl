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
    v = zeros(Int, P)       # Compteur de satellites par plan
    r = rand(1:P, N)
    for k in r
        @inbounds v[k] += 1 # Ajoute un satellite à un plan choisi au hasard
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

"""
    fitness_fixedN(vec, F, i_deg, a, eps_deg; Cmin=75.0, Pbonus=true)

Évalue la qualité d'une constellation candidate.

# Arguments
- vec      : Vecteur de taille P indiquant le nombre de satellites dans chaque plan orbital.
- F        : Paramètre de phasage (Walker-Delta) utilisé par `eval_constellation`.
- i_deg    : Inclinaison orbitale en degrés.
- a        : Demi-grand axe de l'orbite (en mètres), généralement Re + altitude.
- eps_deg  : Angle d'élévation minimal (en degrés) nécessaire pour qu'un satellite couvre un point au sol.

# Paramètres optionnels
- Cmin     : Seuil minimal de couverture acceptable. Si la couverture retournée par `eval_constellation` est < Cmin, une pénalité de -100 est appliquée.
- Pbonus   : Si true, récompense légèrement les vecteurs qui utilisent moins de plans orbitaux (bonus proportionnel au nombre de plans vides).

# Valeur retournée
- fit      : Score de qualité de la constellation.  
             Maximisé lorsque la couverture est haute et que la structure utilise peu de plans.
"""
function fitness_fixedN(vec, F, i_deg, a, eps_deg; Cmin=75.0, Pbonus=true)
    
    key = Tuple(vec)
    if haskey(FITCACHE, key)
        return FITCACHE[key]
    end

    cov, N = eval_constellation(vec, F, i_deg, a, eps_deg; n=10, dlat=6, dlon=6) # Maillage grossier (+ rapide et pas de grande diff dans les valeurs)

    n_used = count(!iszero, vec)
    bonus = Pbonus ? 0.2 * (length(vec) - n_used) : 0 # Si Pbonus=true, on ajoute un bonus en fonction du nombre de plans vides

    fit = cov < Cmin ? cov - 100.0 : cov + bonus # Si cov < Cmin, on applique une forte pénalité
    FITCACHE[key] = fit
    return fit
end

"""
    evolve_vec_fixedN(P, N, F, i_deg, a, eps_deg; popsize=20, generations=30, Cmin=0.0, Pbonus=true)

Algorithme génétique simple pour optimiser la répartition de N satellites sur P plans orbitaux.

# Arguments
- P        : Nombre total de plans orbitaux possibles.
- N        : Nombre total de satellites à répartir sur les P plans.
- F        : Paramètre de phasage (Walker-Delta) utilisé dans `eval_constellation` pour définir la géométrie de la constellation.
- i_deg    : Inclinaison orbitale en degrés.
- a        : Demi-grand axe de l'orbite (en mètres), typiquement Re + altitude.
- eps_deg  : Angle d'élévation minimal (en degrés) pour considérer qu'un satellite couvre un point au sol.

# Paramètres optionnels
- popsize      : Taille de la population de vecteurs candidats (nombre de solutions évaluées par génération).
- generations  : Nombre de générations de l'algorithme génétique (profondeur de la recherche).
- Cmin         : Seuil minimal de couverture ; si la couverture moyenne est < Cmin, une pénalité forte est appliquée dans `fitness`.
- Pbonus       : Si true, ajoute un bonus dans `fitness` pour les configurations qui utilisent moins de plans orbitaux (favorise les plans vides).

# Valeurs de retour
- best_vec : Vecteur de taille P contenant le nombre de satellites par plan pour la meilleure configuration trouvée.
- cov      : Couverture moyenne correspondante (telle que renvoyée par `eval_constellation`, calculée avec une résolution plus fine à la fin).
"""
function evolve_vec_fixedN(P, N, F, i_deg, a, eps_deg; popsize=20, generations=30, Cmin=0.0, Pbonus=true)

    population = [random_vec(P, N) for _ in 1:popsize] # Population initiale
    best_vec = population[1]
    best_fit = fitness_fixedN(best_vec, F, i_deg, a, eps_deg; Cmin=Cmin, Pbonus=Pbonus) # Initialisation de best_fit

    for _ in 1:generations
        fits = [fitness_fixedN(v, F, i_deg, a, eps_deg; Cmin=Cmin, Pbonus=Pbonus) for v in population]
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

    cov, _ = eval_constellation(best_vec, F, i_deg, a, eps_deg; n=100, dlat=1, dlon=1) # Plus fin pour le cov final
    return best_vec, cov
end

"""
    mutate_vec(vec; p_mut=0.3, p_add=0.01)

Applique une mutation sur la configuration orbitale représentée par `vec`.

# Arguments
- `vec  : Vecteur de taille `P` indiquant le nombre de satellites par plan orbital.

# Paramètres optionnels
- `p_mut    : Probabilité de déplacer un satellite d'un plan vers un autre (répété `P` fois).
- `p_add    : Probabilité d'ajouter un satellite sur un plan choisi aléatoirement.

# Retour
- Nouveau vecteur muté, basé sur `vec` (l'original n'est pas modifié).
"""
function mutate_vec(vec; p_mut=0.3, p_add=0.01)
    P = length(vec)
    v = copy(vec)

    for _ in 1:P
        if rand() < p_mut # Proba relativement élevée (on le fait P fois) de modifier la configuration 
            i, j = rand(1:P, 2)
            if v[i] > 0
                v[i] -= 1
                v[j] += 1
            end
        end
    end

    if rand() < p_add # Faible proba d'ajouter un satellite sur un des plans
        k = rand(1:P)
        v[k] += 1
    end

    return v
end

# Dictionnaire pour stocker des configurations et une couverture associé (pour ne pas recalculer ce qui a déjà été calculé)
const FITCACHE = Dict{Tuple{Vararg{Int}}, Float64}()

"""
    fitness(vec, F, i_deg, a, eps_deg; Cmin=95.0, Pbonus=true, Pbonus_coef=0.2, Npenalty=true, Npenalty_coef=0.3)

Fonction de coût qui évalue la qualité d'une configuration en fonction de la couverture et (optionel) du nombre de plans vides et du nombre total de satellites.

# Arguments
- vec         : Vecteur de taille P indiquant le nombre de satellites par plan orbital.
- F           : Paramètre de phasage (Walker-Delta).
- i_deg       : Inclinaison orbitale en degrés.
- a           : Demi-grand axe orbital (en mètres).
- eps_deg     : Angle d'élévation minimal requis pour la visibilité.

# Paramètres optionnels
- Cmin            : Seuil minimal de couverture. Si la couverture retournée est < Cmin, aucune pénalité liée au nombre de satellites n'est prise en compte.
- Pbonus          : Active un bonus proportionnel au nombre de plans orbitaux inutilisés.
- Pbonus_coef     : Coefficient du bonus appliqué par plan vide.
- Npenalty        : Active une pénalité proportionnelle au nombre total de satellites.
- Npenalty_coef   : Coefficient de pénalité par satellite.

# Valeur retournée
- fit : Score de la configuration.
"""
function fitness(vec, F, i_deg, a, eps_deg; Cmin=95.0, Pbonus=true, Pbonus_coef=0.2, Npenalty=true, Npenalty_coef=0.3)

    key = Tuple(vec)
    if haskey(FITCACHE, key)
        return FITCACHE[key]
    end
    
    cov, N = eval_constellation(vec, F, i_deg, a, eps_deg; n=10, dlat=6, dlon=6)

    if cov > Cmin
        # On applique un bonus sur le nombre de plans vides (si Nbonus=true)
        bonus_plans = Pbonus ? Pbonus_coef * (length(vec) - count(!iszero, vec)) : 0
        
        # On applique une pénalité sur N (si Npenalty=true)
        penalty_sat = Npenalty ? Npenalty_coef * N : 0
        
        fit = cov - penalty_sat + bonus_plans
    else
        fit = cov
    end

    FITCACHE[key] = fit
    return fit
end


"""
    evolve_vec(P, N, F, i_deg, a, eps_deg; popsize=20, generations=30, Cmin=0.0, Pbonus=true, Pbonus_coef=0.2, Npenalty=true, Npenalty_coef=0.3, p_mut=0.3, p_add=0.01)

Algorithme génétique optimisant la répartition de `N` satellites sur `P`
plans orbitaux. Chaque solution candidate est un vecteur d'entiers représentant
le nombre de satellites par plan orbital, évalué via la fonction `fitness`.

# Arguments
- P             : Nombre total de plans orbitaux disponibles.
- N             : Nombre total de satellites à répartir sur ces P plans.
- F             : Paramètre de phasage (Walker-Delta) utilisé par `eval_constellation`.
- i_deg         : Inclinaison orbitale en degrés.
- a             : Demi-grand axe de l'orbite (en mètres).
- eps_deg       : Angle d'élévation minimal requis pour qu'un satellite couvre un point au sol.

# Paramètres optionnels
- popsize          : Taille de la population de candidats par génération.
- generations      : Nombre total de générations de l'algorithme génétique.
- Cmin            : Seuil minimal de couverture. Si la couverture retournée est < Cmin, aucune pénalité liée au nombre de satellites n'est prise en compte.
- Pbonus           : Active un bonus dans `fitness` pour les configurations utilisant moins de plans orbitaux.
- Pbonus_coef      : Coefficient du bonus par plan vide.
- Npenalty         : Active une pénalité liée au nombre total de satellites.
- Npenalty_coef    : Coefficient de la pénalité par satellite.
- p_mut            : Probabilité de mutation lors de `mutate_vec`.
- p_add            : Probabilité d'ajouter un satellite lors de `mutate_vec`.

# Valeurs de retour
- best_vec : Vecteur optimal trouvé (taille P), donnant la répartition des satellites.
- cov      : Couverture moyenne finale (résolution fine, n=100, dlat=1, dlon=1).
- N        : Nombre total de satellites correspondant à `best_vec`.
"""
function evolve_vec(P, N, F, i_deg, a, eps_deg; popsize=20, generations=30, Cmin=0.0, Pbonus=true, Pbonus_coef=0.2, Npenalty=true, Npenalty_coef=0.3, p_mut=0.3, p_add=0.01)

    population = [random_vec(P, N) for _ in 1:popsize] # Population initiale
    best_vec = population[1]
    best_fit = fitness(best_vec, F, i_deg, a, eps_deg; Cmin=Cmin, Pbonus=Pbonus, Pbonus_coef=Pbonus_coef, Npenalty=Npenalty, Npenalty_coef=Npenalty_coef)

    for _ in 1:generations
        fits = [fitness(v, F, i_deg, a, eps_deg; Cmin=Cmin, Pbonus=Pbonus, Pbonus_coef=Pbonus_coef, Npenalty=Npenalty, Npenalty_coef=Npenalty_coef) for v in population]
        order = sortperm(fits, rev=true) # Classement par fitness
        elite = population[order[1:clamp(popsize÷4, 1, popsize)]] # On prend que 1/4 des vecteurs (les meilleurs)
        
        if fits[order[1]] > best_fit
            best_fit = fits[order[1]] # Mise à jour de best_fit
            best_vec = elite[1] # Mise à jour de best_vec
        end
        
        newpop = copy(elite)
        while length(newpop) < popsize
            p = elite[rand(1:end)]
            child = mutate_vec(p; p_mut=p_mut, p_add=p_add)   # Mutation d'un vecteur aléatoire
            push!(newpop, child)    # Ajout du mutant
        end
        population = newpop
    end

    cov, N = eval_constellation(best_vec, F, i_deg, a, eps_deg; n=100, dlat=1, dlon=1) # Plus fin pour le cov final
    return best_vec, cov, N
end