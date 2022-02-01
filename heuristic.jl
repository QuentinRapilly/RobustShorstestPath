# greedy_weight algo: optimisation gloutonne avec simple retour arrière sur les pondérations des villes
function greedy_weight(n::Int, s::Int, t::Int, p::Array{Int,1}, S::Int, p_hat::Array{Int,1}, 
    d1::Int, d2::Int, roads::Array{Int64, 2}, d::Vector{Float64}, D::Vector{Float64}, exist_road::Array{Int,2})

    taken_roads = copy(exist_road)

    last_city = -1
    current_city = s
    remaining_cities = Array(zeros(n))
    current_path = Array([s])
    current_path_time = 0

    # permutation permettant de trier la liste des poids
    indices = sortperm(p)

    while current_city != t
        
        # Recherche de la ville avec la ponderation minimale parmi les villes non-explorées
        j=0
        # Tant qu'il reste des villes à étudier, on regarde celle dont le poids est minimum et on l'ajoute au chemin
        while j<n && (taken_roads[current_city][indices[j]] != 1 || sum( p[i]*(1+ p_hat[i]) for i in current_path) > S)
            j = j+1
        end

        if j==n # Pas de ville possible
            deleteat!(current_path, last_city)                                      # On enlève la ville du chemin emprunté
            current_path_time = current_path_time - D[last_city][current_city]      # Decrement du temps mis pour emprunter l'arete consideree
            
            # Réinitialise les routes potentielles partant de current_city
            for j in 1:n
                if taken_roads[current_city][j] == -1
                    taken_roads[current_city][j] = 1
                end
            end
            
            # Traduit dans taken_roads le fait que l'on ne doit plus emprunter cette arete
            taken_roads[last_city][current_city] = -1
            current_city = last_city
            last_city = current_path[-1]
        
        else
            last_city = current_city
            current_city = indices[j]

            # Mise à jour de l'historique des chemins empruntés
            taken_roads[last_city][current_city] = -1
            # Update du temps de parcours
            current_path_time = current_path_time + D[last_city][current_city]
            # Ajout de la ville trouvée au chemin courant
            push!(current_path, current_city)
        end
    end
end
