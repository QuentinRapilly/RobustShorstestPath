using Random
using CPLEX

function init_deltas(n::Int, sparse_D::Array{Float64,2}, d1, d2)
    delta_2 = rand(Float64, n)
    delta_2 = rand()*d2*delta_2/sum(delta_2)
    delta_2 = [min(2,delta_2[i]) for i in 1:n]

    delta_1 = rand(Float64, n, n).*sparse_D
    delta_1 = rand()*d1*delta_1/sum(delta_1)
    delta_1 = [min(sparse_D[i,j],delta_1[i,j]) for i in 1:n, j in 1:n]

    return delta_1, delta_2
end

# greedy_weight algo: optimisation gloutonne avec simple retour arrière sur les pondérations des villes
function greedy_weight(n::Int, s::Int, t::Int, p::Array{Int,1}, S::Int, p_hat::Array{Int,1}, 
    d1::Int, d2::Int, roads::Array{Int64, 2}, d::Array{Float64, 1}, D::Array{Float64, 1}, exist_road::Array{Int,2},
    time_limit::Int, alpha::Float64 = 0.01)

    start = time()

    sparse_d = zeros(n,n)
    sparse_D = zeros(n,n)
    nb_roads = size(d, 1)

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
    delta_1, delta_2 = init_deltas(n, sparse_D, d1, d2)

    while current_city != t && (time() - start < time_limit)
        #println("Taken roads :")
        #display(taken_roads)
        #println("\nCurrent city : ", current_city, " last city : ",last_city)
        #println("Current path :", current_path)
        #readline()

        # Recherche de la ville avec la ponderation minimale parmi les villes non-explorées
        ind_min = 1
        min = Inf
        # Tant qu'il reste des villes à étudier, on regarde celle dont le poids est minimum et on l'ajoute au chemin
        for j in 1:n 
            f_j = (p[j] + delta_2[j]*p_hat[j])+alpha*sparse_d[current_city,j]*(1+delta_1[j])
            if taken_roads[current_city,j]==1 && f_j <min
                min = f_j
                ind_min = j
            end
        end

        #println("Min = ",min," a l'indice :",ind_min)
        if min == Inf || sum(p[i] + delta_2[i]*p_hat[i] for i in current_path) + p[ind_min] + delta_2[ind_min]*p_hat[ind_min] > S
            #println("Pas de chemin")

            if length(current_path)==1 # No possible path.
                # Reinitialise delta_2 pour retester avec une autre valeur
                # println("No possible path.")
                delta_1, delta_2 = init_deltas(n, sparse_D, d1, d2)

                current_city = s
                last_city = s
                taken_roads = copy(exist_road)
                it += 1
                #println("it = ", it)
                
            else
                # Réinitialise les routes potentielles arrivant en current_city
                for i in 1:n
                    if taken_roads[current_city, i] == -1
                        taken_roads[current_city, i] = 1
                    end
                end
                

                # On revient en arriere
                current_city = last_city

                #println("deleted ", current_path[length(current_path)], " from the path")
                deleteat!(current_path, length(current_path)) # On enlève la ville du chemin emprunté
                #println("taille current path : ",length(current_path))
                if length(current_path)==1
                    last_city=s
                else
                    last_city = current_path[length(current_path)-1]
                end
            end
        
        else
            #println("Chemin trouve")
            # Mise à jour de l'historique des chemins empruntés
            taken_roads[current_city, ind_min] = -1
            taken_roads[ind_min, current_city] = -1
            for i in 1:n
                if taken_roads[i, current_city]==1
                    taken_roads[i, current_city] = -1 # On ne passe qu'une seule fois par une ville
                end
            end

            last_city = current_city
            current_city = ind_min
            
            # Ajout de la ville trouvée au chemin courant
            push!(current_path, current_city)
            #println("Added ", current_city, " to the path.")
            #println("current_path = ", current_path)
        end
        it<10000 || break
    end

    if current_city == t
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
        MOI.set(m, MOI.Silent(), true)
        @variable(m, delta_1[1:n, 1:n] >= 0)
        @objective(m, Max, sum(sparse_d[i,j]*(1+delta_1[i,j])*x[i,j] for i in 1:n, j in 1:n))
        @constraint(m, [i in 1:n, j in 1:n], delta_1[i,j] <= sparse_D[i,j])
        @constraint(m, sum(delta_1[i, j] for i in 1:n, j in 1:n) <= d1)
        optimize!(m)

        objective = JuMP.objective_value(m)

        taken_roads = [(i,j) for i in 1:n, j in 1:n if exist_road[i,j]*x[i,j]==1]
        visited_cities = [i for i in 1:n if y[i]==1]

        return taken_roads, visited_cities, objective, 0
    else
        return [], [], -1, 100
    end
end
