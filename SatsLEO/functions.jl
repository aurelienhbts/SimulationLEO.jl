## FONCTIONS dev mais finalement pas si ouf
# Je les garde ici, on sait jamais que ça puisse être utile plus tard

# Va juste augmenter le nombre de satellites et ne va pas trouver de configuration optimale car N va juste augmenter
"""
improve_vec!(vec, F, i_deg, a, eps_deg)

Améliore le vecteur initial afin que la couverture moyenne sur une période augmente. 
"""
function improve_vec!(vec0, F, i_deg, a, eps_deg)
    cov, N = eval_constellation(vec0, F, i_deg, a, eps_deg) # Couverture initiale
    improved = true

    while improved && N < 50 # 50 satellites maximum
        improved = false
        best_gain = 0.0
        best_idx = 0

        # Teste l'ajout d'un satellite dans chaque plan
        for p in 1:length(vec0)
            vec0[p] += 1
            cov2, N2 = eval_constellation(vec0, F, i_deg, a, eps_deg)
            gain = cov2 - cov
            if gain > best_gain
                best_gain = gain
                best_idx = p
            end
            vec0[p] -= 1
        end

        # Ajoute là où le gain est maximal
        if best_idx != 0
            vec0[best_idx] += 1
            cov, N = eval_constellation(vec0, F, i_deg, a, eps_deg)
            improved = true
        end
    end

    return vec0, cov, N
end

# Fonction pas optimale car je dois choisir les poids à la main
"""
biased_vec(P, N; w1, w2, w3)

Construit un vecteur répartissant N satellites entre les P plans orbitaux avec un biais :
si applicable, les trois premiers plans reçoivent plus de satellites (poids w1, w2 & w3),
les autres ont un poids 1.
"""
function biased_vec(P, N; w1=4.0, w2=4.0, w3=4.0)
    weights = ones(Float64, P)
    weights[1] = w1
    if P ≥ 2
        weights[2] = w2
    end
	if P ≥ 3
		weights[3] = w3
	end

    totalw = sum(weights)
    v = floor.(Int, N .* weights ./ totalw)   # Allocation proportionnelle aux poids
    diff = N - sum(v)                         # Correction pour respecter exactement N

    k = 1
    while diff > 0
        v[k] += 1
        diff -= 1
        k = k == P ? 1 : k + 1
    end

    return v
end

# Bien pour initialiser un vecteur mais il faut l'améliorer
"""
balanced_vec(P, N)

Construit un vecteur de longueur P contenant une répartition aussi uniforme que possible
de N satellites entre les P plans orbitaux.
"""
function balanced_vec(P, N)
    base = N ÷ P              # Nombre minimal de satellites par plan
    r = N % P                 # Plans qui recevront un satellite supplémentaire
    v = fill(base, P)         # Répartition uniforme initiale
    for k in 1:r
        v[k] += 1             # Ajout des satellites restants
    end
    return v
end

# Pas mal mais la fonction ne va pas trouver la configuration optimale car ça va dépendre de l'ordre et elle ne passe qu'une seule fois sur chaque couple i,j. 
"""
improve_vec_fixed_N!(vec, F, i_deg, a, eps_deg)

Optimise la répartition de N satellites sur P plans orbitaux (nombre de satellites fixe) en déplaçant un satellite d'un plan vers un autre lorsque cela augmente la couverture moyenne sur une période. L’algorithme répète les déplacements les plus bénéfiques jusqu’à convergence ou jusqu’à 100 itérations.
"""
function improve_vec_fixed_N!(vec, F, i_deg, a, eps_deg)
	
    cov, N = eval_constellation(vec, F, i_deg, a, eps_deg)
    P = length(vec)

    for iter in 1:100
        best_gain = 0.0
        best_i, best_j = 0, 0

        for i in 1:P, j in 1:P
            i == j && continue
            vec[j] == 0 && continue

            vec[i] += 1; vec[j] -= 1 # Test de la config
	            cov2, _ = eval_constellation(vec, F, i_deg, a, eps_deg)
	            gain = cov2 - cov
	            if gain > best_gain # Si ça améliore, on garde i et j pour le refaire après
	                best_gain = gain
	                best_i, best_j = i, j
	            end
            vec[i] -= 1; vec[j] += 1 # On remet comme avant
        end

        if best_gain <= 0 # Si convergence, break
            break
        end

        vec[best_i] += 1; vec[best_j] -= 1 # On refait 
        cov, _ = eval_constellation(vec, F, i_deg, a, eps_deg)
    end

    return vec, cov
end