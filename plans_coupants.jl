using JuMP
using CPLEX

function plans_coupants(n :: Int, s::Int, t::Int, Mat :: Array{Int,2}, p::Array{Int,1}, S::Int, p_hat::Array{Int,1})

    # Create the model
    m = Model(CPLEX.Optimizer)

    # 1 seul thread pour utiliser les callbacks
    MOI.set(m, MOI.NumberOfThreads(), 1)

    ### Variables
    @variable(m, x[1:n, 1:n], Bin)          # vaut 1 si l'on emprunte l'arete ij
    @variable(m, y[1:n], Bin)               # vaut 1 si le trajet passe par la ville d'indice i
    @variable(m, z)
    
    #######################
    ### Problème maître ###
    #######################

    ### Objectif
    @objective(m, Min, z)

    ### Contraintes

    # TODO : rajouter le fait que c'est pour tout delta_1
    @constraint(m, sum(d[i][j]*(1+delta_1[i][j])*x[i][j] for i in 1:n, j in 1:n) <= z)  # Contrainte principale traduisant objectif
    @constraint(m, y[s]==1)                                                             # Chemin passe par s
    @constraint(m, y[t]==1)                                                             # Chemin passe par t
    @constraint(m, [i in 1:n; i!=s], sum(x[j][j] for j in 1:n) == y[i])                 # si on passe par une ville (autre que s), un chemin doit y rentrer
    @constraint(m, [i in 1:n; i!=t], sum(x[i][j] for j in 1:n) == y[i])                 # si on passe par une ville (autre que t), un chemin doit en sortir
    # TODO : trouver le moyen de gérer une contrainte avec Max (via callback ?)
    @constraint(m,  Sum(p[i]*(1+delta_2[i])*y[i] for i in 1:n) <= S)                # somme des poids des villes avec aleas inferieure à un seuil S      

    # U_1_star = {}
    # U_2_star = {}


    #############################
    ### Sous-problème 1 (SPo) ###
    #############################
    function Spo_callback(cb_data::CPLEX.CallbackContext, context_id::Clong)
        ### Objectif
        @objective(m, Max, sum(d[i][j]*(1+delta_1[i][j])*x[i][j] for i in 1:n, j in 1:n))

        ### Contraintes
        @constraint(m, sum(delta_1[i][j] for i  in 1:n, j in 1:n) <= d[1])
        @constraint(m, [i in 1:n, j in 1:n], delta_1[i][j] >= 0) 
        @constraint(m, [i in 1:n, j in 1:n], delta_1[i][j] <= D[i][j])
        
        ### Optimisation
        optimize!(m)
    end


    #############################
    ### Sous-problème 2 (SPj) ###
    #############################
    function Spj_callback(cb_data::CPLEX.CallbackContext, context_id::Clong)
        ### Objectif
        @objective(m, Max, sum(p[i] + delta_2[i]*p_hat[i] for i in 1:n))

        ### Contraintes
        @constraint(m, [i in 1:n], delta_2[i] >= 0)
        @constraint(m, [i in 1:n], delta_2[i] <= 2)

        # Utilisation des callbacks
        MOI.set(m, CPLEX.CallbackFunction(), Spo_callback)
        MOI.set(m, CPLEX.CallbackFunction(), Spj_callback)
    end

    ### Optimisation
    optimize!(m)

end
