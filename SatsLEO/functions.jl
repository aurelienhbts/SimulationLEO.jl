## FONCTIONS dev mais finalement pas si ouf




# Va juste augmenter le nombre de satellites et ne va pas trouver de configuration optimale
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