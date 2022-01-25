using JuMP
using CPLEX

include("plans_coupants.jl")
include("dualisation.jl")
include("reading_files.jl")
#include("heuristic.jl")

PATH_DATA = "data_small/"

function main()

    methods = [plans_coupants, dualisation]#, greedy_weight]

    for file in readdir(PATH_DATA)
        if isfile(PATH_DATA * file) && file != ".directory"
            n, s, t, S, d1, d2, p, ph, Mat, exist_road = read_file(PATH_DATA * file)

            for method in methods
                print("Instance : %s, method : %s\n", file, method)

                start = time()
                x, y = method(n, s, t, p, S, ph, d1, d2, Mat, exist_road)
                t = time() - start
                
                printf("Time taken to proceed instance : %.4f\n", t)

            end

        end
    end
end

main()