using JuMP
using CPLEX

function dualisation(n::Int, s::Int, t::Int, p::Array{Int,1},
    S::Int, p_hat::Array{Int,1}, d1::Int, d2::Int, roads::Array{Int64, 2}, d::Vector{Float64}, D::Vector{Float64}, exist_road::Array{Int,2})

    nb_roads = size(roads,1)
    # Create the model
    m = Model(CPLEX.Optimizer)
    MOI.set(m, MOI.Silent(), true)

    ### Variables
    @variable(m, y[1:n], Bin)
    @variable(m, x[1:n,1:n], Bin)

    @variable(m, alpha_1[1:n,1:n] >= 0)
    @variable(m, alpha_2[1:n] >= 0)
    @variable(m, alpha_1_0 >= 0)
    @variable(m, alpha_2_0 >= 0)

    ### Objective
    @objective(m, Min, sum(x[roads[i,1],roads[i,2]]*d[i] + alpha_1[roads[i,1],roads[i,2]]*D[i] for i in 1:nb_roads) + alpha_1_0*d1)

    ### Constraints
    @constraint(m, [i in 1:nb_roads], alpha_1[roads[i,1], roads[i,2]] + alpha_1_0 >= d[i]*x[roads[i,1], roads[i,2]])
    @constraint(m, [i in 1:n], alpha_2[i] + alpha_2_0 >= p_hat[i]*y[i])
    @constraint(m, sum(p[i]*y[i]+2*alpha_2[i] for i in 1:n) + alpha_2_0*d2 <= S)
    
    @constraint(m, y[s]==1)
    @constraint(m, y[t]==1)


    @constraint(m, [i in 1:n; i!=s], sum(x[j,i]*exist_road[j,i] for j in 1:n) == y[i])
    @constraint(m, [i in 1:n; i!=t], sum(x[i,j]*exist_road[i,j] for j in 1:n) == y[i])
    @constraint(m, [i in 1:n, j in 1:n], x[i,j]+x[j,i] <= 1)

    
    # Optimize the problem
    optimize!(m)

    value_x = JuMP.value.(x)
    taken_roads = [(i,j) for i in 1:n, j in 1:n if exist_road[i,j]*value_x[i,j]==1]

    value_y = JuMP.value.(y)
    visited_cities = [i for i in 1:n if value_y[i]==1]

    obj = JuMP.objective_value(m)

    return taken_roads, visited_cities, obj
end
