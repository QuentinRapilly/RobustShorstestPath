function read_file(file_name :: String)

    
    if isfile(file_name)
        # L’ouvrir
        include(file_name)
        # Lire toutes les lignes d’un fichier
        """data = readlines(myFile) Retourne un tableau de String
        # Pour chaque ligne du fichier
        
        for line in data
            # Afficher la ligne
            println(line)
        end
        # Fermer le fichier
        close(myFile)"""
    end
end