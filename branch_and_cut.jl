using JuMP
using CPLEX

function branch_and_cut(n::Int, s::Int, t::Int, p::Array{Int,1},
    S::Int, p_hat::Array{Int,1}, d1::Int, d2::Int, roads::Array{Int64, 2},
    d::Vector{Float64}, D::Vector{Float64}, exist_road::Array{Int,2}, time_limit::Int)

    nb_roads = size(roads, 1)
    # Create the model
    m = Model(CPLEX.Optimizer)
    set_time_limit_sec(m, time_limit)
    MOI.set(m, MOI.Silent(), true)

    # 1 seul thread pour utiliser les callbacks
    MOI.set(m, MOI.NumberOfThreads(), 1)

    #######################
    ### Problème maître ###
    #######################
    
    ### Variables
    @variable(m, x[1:n, 1:n], Bin)      # vaut 1 si l'on emprunte l'arete ij
    @variable(m, y[1:n], Bin)           # vaut 1 si le trajet passe par la ville d'indice i
    @variable(m, z)

    ### Objectif
    @objective(m, Min, z)

    ### Contraintes
    @constraint(m, sum(d[i]*x[ roads[i,1], roads[i,2] ] for i in 1:nb_roads) <= z)        # Contrainte principale traduisant objectif
    @constraint(m, y[s]==1)                                                 # Chemin passe par s
    @constraint(m, y[t]==1)                                                 # Chemin passe par t
    @constraint(m, [i in 1:n; i!=s], sum(x[j,i]*exist_road[j,i] for j in 1:n) == y[i])     # si on passe par une ville (autre que s), un chemin doit y rentrer
    @constraint(m, [i in 1:n; i!=t], sum(x[i,j]*exist_road[i,j] for j in 1:n) == y[i])     # si on passe par une ville (autre que t), un chemin doit en sortir
    @constraint(m, [i in 1:n, j in 1:n], x[i,j]+x[j,i] <= 1)
    @constraint(m, sum(p[i]*y[i] for i in 1:n) <= S)                        # somme des poids des villes avec aleas inferieure à un seuil S

    function callback(cb_data::CPLEX.CallbackContext, context_id::Clong)

        ### Chargement des variables courantes du problème maître
        if  context_id != CPX_CALLBACKCONTEXT_CANDIDATE
            return 
        end

        CPLEX.load_callback_variable_primal(cb_data, context_id)
        x_val = callback_value.(cb_data, x)
        y_val = callback_value.(cb_data, y)
        z_val = callback_value.(cb_data, z)

        #############################
        ### Sous-problème 1 (SPo) ###
        #############################

        ### Modele pour sous-probleme 1
        m1 = Model(CPLEX.Optimizer)
        MOI.set(m1, MOI.Silent(), true)

        # Variable
        @variable(m1, delta_1[1:n, 1:n] >= 0)    # aleas sur le temps de trajet

        ### Objectif
        @objective(m1, Max, sum(d[i]*(1+delta_1[ roads[i,1], roads[i,2] ])*x_val[roads[i,1], roads[i,2]] for i in 1:nb_roads))

        ### Contraintes
        @constraint(m1, sum(delta_1[roads[i,1], roads[i,2]] for i in 1:nb_roads) <= d1)
        @constraint(m1, [i in 1:nb_roads], delta_1[roads[i,1], roads[i,2]] <= D[i])
        
        ### Optimisation
        optimize!(m1)

        ### Valeur optimale de delta_1 retournée par le solveur
        delta_1_star = JuMP.value.(delta_1)
        tmp1 = sum(d[i]*(1+delta_1_star[roads[i,1], roads[i,2]])*x_val[roads[i,1], roads[i,2]] for i in 1:nb_roads)
        if  tmp1 > z_val
            # Ajout de la contrainte au problème maître
            cstr = @build_constraint( sum(d[i]*(1+delta_1_star[roads[i,1], roads[i,2]])*x[roads[i,1], roads[i,2]] for i in 1:nb_roads) <= z)
            MOI.submit(m, MOI.LazyConstraint(cb_data), cstr)
        end

        #############################
        ### Sous-problème 2 (SPj) ###
        #############################

        ### Modele pour le sous-probleme 2
        m2 = Model(CPLEX.Optimizer)
        MOI.set(m2, MOI.Silent(), true)

        # Variable
        @variable(m2, delta_2[1:n] >= 0)       # aleas sur la ponderation des villes

        ### Objectif
        @objective(m2, Max, sum((p[i] + delta_2[i]*p_hat[i])*y_val[i] for i in 1:n))

        ### Contraintes
        @constraint(m2, [i in 1:n], delta_2[i] <= 2)
        @constraint(m2, sum(delta_2[i] for i in 1:n) <= d2)        # limite sur somme aleas ponderation des villes

        ### Optimisation
        optimize!(m2)

        ### Valeur optimale de delta_1 retournée par le solveur
        delta_2_star = JuMP.value.(delta_2)

        tmp2 = sum((p[i]+delta_2_star[i]*p_hat[i])*y_val[i] for i in 1:n)
        if  tmp2 > S
            # Ajout de la contrainte au problème maître
            cstr = @build_constraint( sum((p[i]+delta_2_star[i]*p_hat[i])*y[i] for i in 1:n) <= S)
            MOI.submit(m, MOI.LazyConstraint(cb_data), cstr)
        end
    end


    MOI.set(m, CPLEX.CallbackFunction(), callback)

    ### Optimisation
    optimize!(m)

    value_x = JuMP.value.(x)
    taken_roads = [(i,j) for i in 1:n, j in 1:n if exist_road[i,j]*value_x[i,j]==1]

    value_y = JuMP.value.(y)
    visited_cities = [i for i in 1:n if value_y[i]==1]
    obj = JuMP.objective_value(m)

    GAP = MOI.get(m, MOI.RelativeGap())

    return taken_roads, visited_cities, obj, GAP
end
