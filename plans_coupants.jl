using JuMP
using CPLEX
using Flux



function plans_coupants(n::Int, s::Int, t::Int, p::Array{Int,1},
    S::Int, p_hat::Array{Int,1}, d1::Int, d2::Int, roads::Array{Int64, 2},
    d::Vector{Float64}, D::Vector{Float64}, exist_road::Array{Int,2})

    list_delta_1 = zeros(Float64, 1, n, n)
    list_delta_2 = zeros(Float64,1,n)

    nb_roads = size(roads, 1)

    x=-1
    y=-1
    z=-1

    while true
        added_cstr = false
        x, y, z = master_problem(n, s, t, p, S, p_hat, roads, d, D, exist_road, list_delta_1, list_delta_2)
        delta_1 = SPo(x,  n, roads, d, D, d1, exist_road)
        delta_2 = SPj(y, n, p, p_hat, d2, roads, exist_road)

        if sum(d[i]*(1+delta_1[roads[i,1], roads[i,2]])*x[roads[i,1], roads[i,2]] for i in 1:nb_roads) > z
            list_delta_1 = vcat(list_delta_1,Flux.unsqueeze(delta_1,1))
            added_cstr = true
        end

        if sum((p[i]+delta_2[i]*p_hat[i])*y[i] for i in 1:n) > S
            list_delta_2 = vcat(list_delta_2, Flux.unsqueeze(delta_2,1))
            added_cstr = true
        end

        (! added_cstr) && break

    end

    taken_roads = [(i,j) for i in 1:n, j in 1:n if exist_road[i,j]*x[i,j]==1]
    visited_cities = [i for i in 1:n if y[i]==1]

    return taken_roads, visited_cities, z

end

function master_problem(n::Int, s::Int, t::Int, p::Array{Int,1},
    S::Int, p_hat::Array{Int,1}, roads::Array{Int, 2},
    d::Vector{Float64}, D::Vector{Float64}, exist_road::Array{Int,2},
    list_delta_1::Array{Float64,3}, list_delta_2::Array{Float64,2})

    n_delta_1 = size(list_delta_1, 1)
    n_delta_2 = size(list_delta_2, 1)

    nb_roads = size(roads, 1)
    # Create the model
    m = Model(CPLEX.Optimizer)
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
    @constraint(m, [k in 1:n_delta_1], sum(d[i]*(1+list_delta_1[k,roads[i,1], roads[i,2]])*x[roads[i,1], roads[i,2]] for i in 1:nb_roads) <= z)        # Contrainte principale traduisant objectif
    @constraint(m, y[s]==1)                                                 # Chemin passe par s
    @constraint(m, y[t]==1)                                                 # Chemin passe par t
    @constraint(m, [i in 1:n; i!=s], sum(x[j,i]*exist_road[j,i] for j in 1:n) == y[i])     # si on passe par une ville (autre que s), un chemin doit y rentrer
    @constraint(m, [i in 1:n; i!=t], sum(x[i,j]*exist_road[i,j] for j in 1:n) == y[i])     # si on passe par une ville (autre que t), un chemin doit en sortir
    @constraint(m, [i in 1:n, j in 1:n], x[i,j]+x[j,i] <= 1)
    @constraint(m, [k in 1:n_delta_2],sum((p[i]+list_delta_2[k,i]*p_hat[i])*y[i] for i in 1:n) <= S)

    ### Optimisation
    optimize!(m)
    return JuMP.value.(x), JuMP.value.(y), JuMP.value.(z)
end

function SPo(x::Array{Float64,2}, n::Int, roads::Array{Int64, 2},
    d::Vector{Float64}, D::Vector{Float64}, d1::Int, exist_road::Array{Int,2})

    ### Modele pour sous-probleme 1
    m1 = Model(CPLEX.Optimizer)
    MOI.set(m1, MOI.Silent(), true)

    nb_roads = size(roads, 1)

    # Variable
    @variable(m1, delta_1[1:n, 1:n] >= 0)    # aleas sur le temps de trajet

    ### Objectif
    @objective(m1, Max, sum(d[i]*(1+delta_1[ roads[i,1], roads[i,2] ])*x[roads[i,1], roads[i,2]] for i in 1:nb_roads))

    ### Contraintes
    @constraint(m1, sum(delta_1[roads[i,1], roads[i,2]] for i in 1:nb_roads) <= d1)
    @constraint(m1, [i in 1:nb_roads], delta_1[roads[i,1], roads[i,2]] <= D[i])
    
    ### Optimisation
    optimize!(m1)
    delta_1_star = JuMP.value.(delta_1)

    return delta_1_star

end

function SPj(y::Array{Float64,1}, n::Int, p::Array{Int,1}, p_hat::Array{Int,1}, 
    d2::Int, roads::Array{Int64, 2}, exist_road::Array{Int,2})

    ### Modele pour le sous-probleme 2
    m2 = Model(CPLEX.Optimizer)
    MOI.set(m2, MOI.Silent(), true)

    nb_roads = size(roads, 1)

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

    return delta_2_star

end