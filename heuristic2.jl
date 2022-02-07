using Random
using CPLEX

function greedy_weight(n::Int, s::Int, t::Int, p::Array{Int,1}, S::Int, p_hat::Array{Int,1}, 
    d1::Int, d2::Int, roads::Array{Int64, 2}, d::Array{Float64, 1}, D::Array{Float64, 1}, exist_road::Array{Int,2},
    time_limit::Int, alpha::Float64 = 0.001)

    start = time()

    sparse_d = zeros(n,n)
    nb_roads = size(d, 1)

    for i in 1:nb_roads
        sparse_d[roads[i,1],roads[i,2]] = d[i]
    end

    mult_factor = 1

    max_ph = -Inf
    for tmp in p_hat
        if tmp>max_ph
            max_ph = tmp
        end
    end

    taken_roads = copy(exist_road)
    for i in 1:n 
        taken_roads[i,s]=0
    end

    it = 0
    last_city = -1
    current_city = s
    current_path = Array([s])


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
        j = 1
        while j<n 
            j != t || taken_roads[current_city,j]!=1 || break
            f_j = p[j] + alpha*sparse_d[current_city,j]
            if taken_roads[current_city,j]==1 && f_j <min
                min = f_j
                ind_min = j
            end
            j += 1
        end

        if j<n
            ind_min = j
        end

        #println("Min = ",min," a l'indice :",ind_min)
        sum_path = sum(p[i] for i in current_path) + p[ind_min]
        reduced_S = S - d2*max_ph*mult_factor
        if min == Inf || sum_path > reduced_S
            #println("Pas de chemin")

            if length(current_path)==1 # No possible path.
                # Reinitialise delta_2 pour retester avec une autre valeur
                #println("No possible path.")
                mult_factor = mult_factor*0.99

                current_city = s
                last_city = s
                taken_roads = copy(exist_road)
        
                # println("Valeur de reduced S : ",reduced_S)
                # readline()
                for i in 1:n 
                    taken_roads[i,s]=0
                end
                it += 1
                # print("it = ", it, "\r")
                
            else
                # Réinitialise les routes potentielles arrivant en current_city
                for i in 1:n
                    if taken_roads[i, current_city] == -1
                        taken_roads[i, current_city] = 1
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

            last_city = current_city
            current_city = ind_min

            taken_roads[last_city, current_city] = -2
            taken_roads[current_city, last_city] = -1

            for i in 1:n
                if taken_roads[i, current_city] == 1
                    taken_roads[i, current_city] = -1 # On ne passe qu'une seule fois par une ville
                end
            end

            
            # Ajout de la ville trouvée au chemin courant
            push!(current_path, current_city)
            #println("Added ", current_city, " to the path.")
            #println("current_path = ", current_path)
        end
        it<500 || break
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

        m2 = Model(CPLEX.Optimizer)
        MOI.set(m2, MOI.Silent(), true)

        # Variable
        @variable(m2, delta_2[1:n] >= 0)       # aleas sur la ponderation des villes

        ### Objectif
        @objective(m2, Max, sum((p[i] + delta_2[i]*p_hat[i])*y[i] for i in 1:n))

        ### Contraintes
        @constraint(m2, [i in 1:n], delta_2[i] <= 2)
        @constraint(m2, sum(delta_2[i] for i in 1:n) <= d2)        # limite sur somme aleas ponderation des villes

        ### Optimisation
        optimize!(m2)

        ### Valeur optimale de delta_1 retournée par le solveur
        delta_2_star = JuMP.value.(delta_2)

        tmp2 = sum((p[i]+delta_2_star[i]*p_hat[i])*y[i] for i in 1:n)
        if  tmp2 > S
            return [], [], -1, 100
        else 

            m = Model(CPLEX.Optimizer)
            MOI.set(m, MOI.Silent(), true)
            @variable(m, delta_1[1:n, 1:n] >= 0)
            @objective(m, Max, sum(d[i]*(1+delta_1[roads[i,1],roads[i,2]])*x[roads[i,1],roads[i,2]] for i in 1:nb_roads))
            @constraint(m, [i in 1:nb_roads], delta_1[roads[i,1],roads[i,2]] <= D[i])
            @constraint(m, sum(delta_1[roads[i,1], roads[i,2]] for i in 1:nb_roads) <= d1)
            optimize!(m)

            objective = JuMP.objective_value(m)

            taken_roads = [(i,j) for i in 1:n, j in 1:n if exist_road[i,j]*x[i,j]==1]
            visited_cities = [i for i in 1:n if y[i]==1]

            return taken_roads, visited_cities, objective, 0
        end

    else
        return [], [], -1, 100
    end
end
