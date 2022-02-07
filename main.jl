using JuMP
using CPLEX

include("branch_and_cut.jl")
include("dualisation.jl")
include("reading_files.jl")
include("heuristic.jl")
include("plans_coupants.jl")

PATH_DATA = "data/"
PATH_RES = "res/"

function main(start_at_file)

    methods = [plans_coupants, dualisation, branch_and_cut]# ,greedy_weight]

    time_limit = 180

    csv_file = open(PATH_RES*"output.csv", "w")
    println(csv_file, "instance,time_PC,obj_PC,gap_PC,time_dual,obj_dual,gap_dual,time_branch,obj_branch,gap_branch")
    close(csv_file)

    files_array = sort_instances_by_size(PATH_DATA)

    for file in files_array[start_at_file:size(files_array,1)]
        if isfile(PATH_DATA * file) && file != ".directory"
            n, s, t, S, d1, d2, p, ph, sparse_roads, sparse_d, sparse_D, exist_road = read_line_by_line(PATH_DATA * file)

            instance = split(file,".")
            p1 = instance[1]
            instance = split(instance[1],"_")[1]*"-"*instance[2]
            output_file = open(PATH_RES*instance*".res","w")

            csv_line = instance*","

            for method in methods
                println("Instance : ",file, ", method : ", method)

                println(output_file, "Method : ", method)

                start = time()
                x, y, obj, gap = method(n, s, t, p, S, ph, d1, d2, sparse_roads, sparse_d, sparse_D, exist_road, time_limit)   
                elapsed_time = time() - start
                
                println("Time taken to proceed instance : ", elapsed_time)
                println(output_file, "Time : ", elapsed_time)
                println(output_file, "x = ", x)
                println(output_file, "y = ", y)
                println(output_file, "obj = ", obj)
                println(output_file,' ')

                csv_line = csv_line*string(round(elapsed_time,digits=2))*","*string(floor(obj))*","*string(round(gap,digits = 1))*","

            end

            csv_file = open(PATH_RES*"output.csv", "a")
            println(csv_file, strip.(csv_line,','))
            close(csv_file)

            close(output_file)
        end
    end

end

main(1)