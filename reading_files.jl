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