using JuMP
using CPLEX

include("branch_and_cut.jl")
include("dualisation.jl")
include("reading_files.jl")
include("heuristic.jl")
include("plans_coupants.jl")

PATH_DATA = "data_small/"
PATH_RES = "res/"

function main()

    methods = [plans_coupants, dualisation, branch_and_cut]#, greedy_weight]

    for file in readdir(PATH_DATA)
        if isfile(PATH_DATA * file) && file != ".directory"
            n, s, t, S, d1, d2, p, ph, sparse_roads, sparse_d, sparse_D, exist_road = read_line_by_line(PATH_DATA * file)
            
            instance = split(file,".")
            instance = instance[1]*"_"*instance[2]
            output_file = open(PATH_RES*instance*".res","w")

            for method in methods
                println("Instance : ",file, ", method : ", method)

                println(output_file, "Method : ", method)

                start = time()
                x, y, obj = method(n, s, t, p, S, ph, d1, d2, sparse_roads, sparse_d, sparse_D, exist_road)
                elapsed_time = time() - start
                
                println("Time taken to proceed instance : ", elapsed_time)

                println(output_file, "Time : ", elapsed_time)
                println(output_file, "x = ",x)
                println(output_file, "y = ",y)
                println(output_file, "obj = ", obj)
                println(output_file,' ')

            end
            close(output_file)
        end
    end
end

main()