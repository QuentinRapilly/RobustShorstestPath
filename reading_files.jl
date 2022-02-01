function read_file(file_name :: String)
    if isfile(file_name)
        # Ouverture du fichier
        include(file_name)
        exist_road = zeros(Int,n,n)
        nb_roads = size(Mat,1)

        # Récupération des routes de l'instance (conversion en int pour indexation car matrice de float)
        sparse_roads = zeros(Int, nb_roads, 2)
        for i in 1:nb_roads
            sparse_roads[i,1] = convert(Int,Mat[i,1]) 
            sparse_roads[i,2] = convert(Int,Mat[i,2])
        end

        # Récupération des d
        sparse_d = [Mat[i, 3] for i in 1:nb_roads]
        # Récupération des D
        sparse_D = [Mat[i, 4] for i in 1:nb_roads]

        for i in 1:nb_roads
            exist_road[sparse_roads[i,1], sparse_roads[i,2]] = 1
        end
        return(n,s,t,S,d1,d2,p,ph,sparse_roads,sparse_d,sparse_D,exist_road)
    end
end

function read_line_by_line(file_name::String)

    if isfile(file_name)

        f = open(file_name)

        # Parsing du début de fichier
        n = parse(Int, last(split(readline(f), " ")))
        s = parse(Int, last(split(readline(f), " ")))
        t = parse(Int, last(split(readline(f), " ")))
        S = parse(Int, last(split(readline(f), " ")))
        d1 = parse(Int, last(split(readline(f), " ")))
        d2 = parse(Int, last(split(readline(f), " ")))
        line = readline(f)
        p = split(line[6:length(line)-1], ", ")
        p = [parse(Int, elem) for elem in p]
        line = readline(f)
        p_hat = split(line[7:length(line)-1], ", ")
        p_hat = [parse(Int, elem) for elem in p_hat]

        # Parsing de la matrice Mat
        sparse_roads = []
        sparse_d = Array{Float64, 1}([])
        sparse_D = Array{Float64, 1}([])
        readline(f) # > "Mat = ["
        while !eof(f)
            line = readline(f)
            arr = split(line[1:length(line)-1], " ") # Removing the ";"
            push!(sparse_roads, [parse(Int, arr[1]) parse(Int, arr[2])])
            push!(sparse_d, parse(Float64, arr[3]))
            push!(sparse_D, parse(Float64, arr[4]))
        end

        roads = zeros(Int, size(sparse_roads, 1), 2)
        for i in 1:size(sparse_roads, 1)
            roads[i, 1:2] = sparse_roads[i]
        end

        # Remplissage de exist_road
        nb_roads = size(sparse_roads, 1)
        exist_road = zeros(Int, n, n)
        for i in 1:nb_roads
            exist_road[ roads[i,1], roads[i,2] ] = 1
        end

        return(n,s,t,S,d1,d2,p,p_hat,roads,sparse_d,sparse_D,exist_road)
    end
end
