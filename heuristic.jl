using Random
using CPLEX

# greedy_weight algo: optimisation gloutonne avec simple retour arrière sur les pondérations des villes
function greedy_weight(n::Int, s::Int, t::Int, p::Array{Int,1}, S::Int, p_hat::Array{Int,1}, 
    d1::Int, d2::Int, roads::Array{Int64, 2}, d::Array{Float64, 1}, D::Array{Float64, 1}, exist_road::Array{Int,2}, time_limit::Int)

    sparse_d = zeros(n,n)
    sparse_D = zeros(n,n)
    nb_roads = length(d)
    for i in 1:nb_roads
        sparse_D[roads[i,1],roads[i,2]] = D[i]
        sparse_d[roads[i,1],roads[i,2]] = d[i]
    end

    taken_roads = copy(exist_road)

    it = 0
    last_city = -1
    current_city = s
    current_path = Array([s])

    # Pire des cas sur l'alea de la ponderation des villes
    delta_2 = rand(Float64, n)
    delta_2 = d2*delta_2/sum(delta_2)

    # permutation permettant de trier la liste des poids
    indices = sortperm(p + delta_2.*p_hat)
    # println("indices = ", indices)

    while current_city != t && (time() - start < time_limit)

        # Recherche de la ville avec la ponderation minimale parmi les villes non-explorées
        j=1
        # Tant qu'il reste des villes à étudier, on regarde celle dont le poids est minimum et on l'ajoute au chemin
        while j<n && taken_roads[current_city, indices[j]] != 1
            # println("taken_roads : ", taken_roads[current_city, indices[j]])
            # println("somme pondération : ", sum( p[i]*(1+ delta_2[i]*p_hat[i] ) for i in current_path))
            j = j+1
        end

        if j == n || sum( p[i] + delta_2[i]*p_hat[i] for i in current_path) + p[indices[j]] + delta_2[indices[j]]*p_hat[indices[j]] > S
            # println("No further path possible.")

            if current_path==[s] # No possible path.
                # Reinitialise delta_2 pour retester avec une autre valeur
                # println("No possible path.")
                delta_2 = rand(Float64, n)
                delta_2 = d2*delta_2/sum(delta_2)
                current_city = s
                last_city = -1
                taken_roads = copy(exist_road)
                it += 1
                println("it = ", it)
                if it > 10000
                    println("Pas de solution trouvée.")
                    return [-1], [-1], -1
                end

            else
                # Réinitialise les routes potentielles arrivant en current_city
                for i in 1:n
                    if taken_roads[i, current_city] == -1
                        taken_roads[i, current_city] = 1
                    end
                end

                # Traduit dans taken_roads le fait que l'on ne doit plus emprunter cette arete
                taken_roads[last_city, current_city] = -1

                # On revient en arriere
                current_city = last_city

                println("deleted ", current_path[length(current_path)], " from the path")
                deleteat!(current_path, length(current_path)) # On enlève la ville du chemin emprunté
                last_city = current_path[length(current_path)]

                println("current_path = ", current_path) 
            end
        
        else
            # Mise à jour de l'historique des chemins empruntés
            taken_roads[current_city, indices[j]] = -1
            for i in 1:n
                if i!=current_city && taken_roads[i, current_city]==1
                    taken_roads[i, current_city] = -1 # On ne passe qu'une seule fois par une ville
                end
            end

            last_city = current_city
            current_city = indices[j]
            
            # Ajout de la ville trouvée au chemin courant
            push!(current_path, current_city)
            println("Added ", current_city, " to the path.")
            println("current_path = ", current_path)
        end
    end

    # Reconstruction de la solution y
    y = zeros(n)
    for i in 1:length(current_path)
        for j in 1:n
            if current_path[i] == j
                y[j] = 1
            end
        end
    end

    # Reconstruction de x à partir de y
    x = zeros(n, n)
    for i in 1:n
        for j in 1:n
            if y[i]==1 && y[j]==1 && i!=j
                x[i,j] = 1
            end
        end
    end

    m = Model(CPLEX.Optimizer)
    @variable(m, delta_1[1:n, 1:n] >= 0)
    @objective(m, Max, sum(sparse_d[i,j]*(1+delta_1[i,j])*x[i,j] for i in 1:n, j in 1:n))
    @constraint(m, [i in 1:n, j in 1:n], delta_1[i,j] <= sparse_D[i,j])
    @constraint(m, sum(delta_1[i, j] for i in 1:n, j in 1:n) <= d1)
    optimize!(m)

    objective = JuMP.objective_value(m)

    taken_roads = [(i,j) for i in 1:n, j in 1:n if exist_road[i,j]*x[i,j]==1]
    visited_cities = [i for i in 1:n if y[i]==1]

    return taken_roads, visited_cities, objective
end
