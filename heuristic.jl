# greedy_weight algo:
#### 1. Initialisation : 
####    - On part de s (Start), on charge les villes adjacentes à s dans G
####    - On regarde la ville qui coûte le moins cher en terme de poids (min_{i \in V(s)} p_i) 
####      où V(s) sont les i tq la ville i est voisine à s
####    - On vérifie que la ville regardée satisfait les contraintes du problème (i.e. ne viole pas la contrainte sur la
####      pondération des villes, et aussi sur la longueur du trajet), on prend la première qui satisfait ces contraintes
####
#### 2. Algorithme
####    - Pour chaque ville courante v, on charge les villes adjacentes à v dans G
####    - On regarde la ville qui coûte le moins cher en terme de poids (min_{i \in V(s)} p_i)
####    - On vérifie que la ville regardée satisfait les contraintes du problème (i.e. ne viole pas la contrainte sur la
####      pondération des villes, et aussi sur la longueur du trajet), on prend la première qui satisfait ces contraintes
####     (simple retour arrière si contraintes pas satisfaites)

function greedy_weight(s::Int, t::Int, G::Array{Int,2}, d::Array{Int,2}, D::Array{Int,2}, p::Array{Int,1}, S::Int, p_hat::Array{Int,1})

    n = size(G, 1)
    last_city = s
    current_city = s
    remaining_cities = Array([i for i in 1:n if i!=s])
    current_path = Array([s]) 
    current_path_time = 0

    # permutation permettant de trier la liste des poids
    indices = sortperm(p)

    while current_city != t

        # Recherche de la ville avec la ponderation minimale parmi les villes non-explorées
        ind_min = indices[0]

        i=0
        remaining_cities_current = Array([i for i in 1:n if !in(current_path, i)])
        while i < n && i in remaining_cities_current
            if G[current_city][indices[i]] == 1
                last_city = current_city
                ind_min = indices[i]
                pop!(remaining_cities_current, current_city)
                break
            else
                i = i+1
            end
        end

        if i==n || isempty(remaining_cities_current) # Pas de ville restante
            deleteat!(current_path, last_city)
            current_path_time = current_path_time - D[last_city][current_city]
            current_city = last_city
        end

        # # Verification des contraintes pour le chemin comportant la ville candidate (ind_min)
        # #TODO : ajouter les contraintes à vérifier pour que la solution en cours soit faisable pour le problème
        # if contrainte_ok 
        #     # Update du temps et de la ville courante
        #     current_path_time = current_path_time + D[current_city][indices[i]]
        #     current_city = ind_min
        #     # Ajout de la ville trouvée au chemin courant
        #     push!(current_path, current_city)
        #     # Enleve la ville visitée et ajoutee au chemin de la liste des villes restantes
        #     deleteat!(remaining_cities, ind_min)
        # end

        # Update du temps et de la ville courante
        current_path_time = current_path_time + D[current_city][indices[i]]
        current_city = ind_min
        # Ajout de la ville trouvée au chemin courant
        push!(current_path, current_city)

    end
end
