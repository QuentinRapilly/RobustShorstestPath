function read_file(file_name::String)
    
    if isfile(file_name)
        # L’ouvrir
        print(file_name)
        include(file_name)
        exist_road = zeros(Int,n,n)
        for r in Mat
            exist_road[r[0],r[1]] = 1
        end
        return(n,s,t,S,d1,d2,p,ph,Mat,exist_road)
    end
end