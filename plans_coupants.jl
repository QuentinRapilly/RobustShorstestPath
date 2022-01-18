using JuMP
using CPLEX

function plans_coupants(n :: Int, s::Int, t::Int, Mat :: Array{Int,2}, p::Array{Int,1}, S::Int, p_hat::Array{Int,1})

    # Create the model
    m = Model(CPLEX.Optimizer)

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
    @constraint(m, sum(d[i][j]*x[i][j] for i in 1:n, j in 1:n) <= z)        # Contrainte principale traduisant objectif
    @constraint(m, y[s]==1)                                                 # Chemin passe par s
    @constraint(m, y[t]==1)                                                 # Chemin passe par t
    @constraint(m, [i in 1:n; i!=s], sum(x[j][j] for j in 1:n) == y[i])     # si on passe par une ville (autre que s), un chemin doit y rentrer
    @constraint(m, [i in 1:n; i!=t], sum(x[i][j] for j in 1:n) == y[i])     # si on passe par une ville (autre que t), un chemin doit en sortir
    @constraint(m, Sum(p[i]*y[i] for i in 1:n) <= S)                        # somme des poids des villes avec aleas inferieure à un seuil S


    #############################
    ### Sous-problème 1 (SPo) ###
    #############################
    function Spo_callback(cb_data::CPLEX.CallbackContext, context_id::Clong)

        ### Modele pour sous-probleme 1
        m1 = Model(CPLEX.Optimizer)

        # Variable
        @variable(m1, delta_1[1:n, 1:n] >= 0)    # aléas sur le temps de trajet

        ### Objectif
        @objective(m1, Max, sum(d[i][j]*(1+delta_1[i][j])*x[i][j] for i in 1:n, j in 1:n))

        ### Contraintes
        @constraint(m1, sum(delta_1[i][j] for i  in 1:n, j in 1:n) <= d[1])
        @constraint(m1, [i in 1:n, j in 1:n], delta_1[i][j] <= D[i][j])
        
        ### Optimisation
        optimize!(m1)

        ### Chargement des variables courantes du problème maître
        CPLEX.load_callback_variable_primal(cb_data, context_id)
        z_val = callback_value(cb_data, z)
        x_val = callback_value(cb_data, x)

        ### Valeur optimale de delta_1 retournée par le solveur
        delta_1_star = CPLEX.get_solution(m1)

        if sum(d[i][j]*(1+delta_1_star[i][j])*x_val[i][j] for i in 1:n, j in 1:n) > z_val
            # Ajout de la contrainte au problème maître
            cstr = @build_constraint( sum(d[i][j]*(1+delta_1_star[i][j])*x[i][j] for i in 1:n, j in 1:n) <= z)
            MOI.submit(m, MOI.UserCut(cb_data), cstr)
        end
    end


    #############################
    ### Sous-problème 2 (SPj) ###
    #############################
    function Spj_callback(cb_data::CPLEX.CallbackContext, context_id::Clong)
        ### Modele pour le sous-probleme 2
        m2 = Model(CPLEX.Optimizer)

        # Variable
        @variable(m2, delta_2[1:n] >= 0)                             # aleas sur la ponderation des villes

        ### Objectif
        @objective(m2, Max, sum(p[i] + delta_2[i]*p_hat[i] for i in 1:n))

        ### Contraintes
        @constraint(m2, [i in 1:n], delta_2[i] <= 2)
        @constraint(m2, sum(delta_2[i] for i in 1:n) <= d[2])        # limite sur somme aleas ponderation des villes

        ### Optimisation
        optimize!(m2)

        ### Chargement des variables courantes du problème maître
        CPLEX.load_callback_variable_primal(cb_data, context_id)
        y_val = callback_value(cb_data, y)

        ### Valeur optimale de delta_1 retournée par le solveur
        delta_2_star = CPLEX.get_solution(m2)

        if Sum((p[i]+delta_2_star[i]*p_hat[i])*y_val[i] for i in 1:n) > S
            # Ajout de la contrainte au problème maître
            cstr = @build_constraint( Sum((p[i]+delta_2_star[i]*p_hat[i])*y[i] for i in 1:n) <= S)
            MOI.submit(m, MOI.UserCut(cb_data), cstr)
        end
    end


    MOI.set(m, CPLEX.CallbackFunction(), Spo_callback)
    MOI.set(m, CPLEX.CallbackFunction(), Spj_callback)

    ### Optimisation
    optimize!(m)

end
