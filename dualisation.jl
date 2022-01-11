using JuMP
using CPLEX

function dualisation(s::Int, t::Int, G::Array{Int,2}, d::Array{Int,2}, D::Array{Int,2}, p::Array{Int,1},
    S::Int, p_hat::Array{Int,1}, d1::Int, d2::Int)

    n = size(G, 1)

    # Create the model
    m = Model(CPLEX.Optimizer)

    ### Variables
    @variable(m, y[1:n], Bin)
    @variable(m, x[1:n,1:n], Bin)

    @variable(m, alpha_1[1:n,1:n] >= 0)
    @variable(m, alpha_2[1:n] >= 0)
    @variable(m, alpha_1_0 >= 0)
    @variable(m, alpha_2_0 >= 0)

    ### Objective
    @objective(m, Min, sum(x[i,j]*d[i,j] + alpha_1[i,j]*D[i,j] for i in 1:n, j in 1:n) + alpha_1_0*d1)

    ### Constraints
    @constraint(m, [i in 1:n, j in 1:n], alpha_1[i,j] + alpha_1_0 >= d[i,j]*x[i,j])
    @constraint(m, [i in 1:n], alpha_2[i] + alpha_2_0 >= p[i]*y[i])
    @constraint(m, sum(p[i]*y[i]+2*alpha_2[i]*p_hat[i] for i in 1:n) + alpha_2_0*d2 <= S)
    
    @constraint(m, y[s]==1)
    @constraint(m, y[t]==1)

    for i in 1:n 
        if i!=s
            @constraint(m, sum(x[j,i] for j in 1:n) == y[i])
        end
        if i!=t 
            @constraint(m, sum(x[i,j] for j in 1:n) == y[i])
        end 
    end
    # Optimize the problem
    optimize!(m)



end