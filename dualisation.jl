using JuMP
using CPLEX

function dualisation(s::Int, t::Int, G::Array{Int,2}, d::Array{Int,2}, D::Array{Int,2}, p::Array{Int,1}, S::Int, p_hat::Array{Int,1})

    n = size(G, 1)

    # Create the model
    m = Model(CPLEX.Optimizer)

    ### Objective
    @objective(m, Max, 1)

    ### Variables
    @variable(m, x[1:n, 1:n, 1:n], Bin)

    ### Constraints

    # Optimize the problem
    optimize!(m)



end