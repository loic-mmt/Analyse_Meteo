using NCDatasets # Pour le chargement des datasets .nc
using Glob #Pour pouvoir chercher des noms de fichiers dans un répèrtoire facillement.
using Statistics #Pour utiliser moyenne (mean) sur des matrices/dataset
using StatsPlots
using Dates # Pour utiliser le format "date" à l'intérieur des fichiers .nc
using GLM # Pour analyses de tendance (modèle linéaire simple)
using DataFrames 
using Base.Threads 
using RollingFunctions # Pour la création simple et rapide de vecteurs pondérés -> courbes pondérés

# Assignation des différents chemins utiles
data_folder_precise = "data/raw_monthly_combined/precise"
data_folder_basic ="data/raw_monthly_combined/basic"
weight_prop_basic_file = "data/masks/weights_prop_basic.nc"
weight_prop_precise_file = "data/masks/weights_prop_precise.nc"
weight_bool_basic_file = "data/masks/weights_bool_basic.nc"
weight_bool_precise_file = "data/masks/weights_bool_precise.nc"

# Importation des poids
ds_p_b = NCDataset(weight_prop_basic_file)
# ds_p_p = NCDataset(weight_prop_precise_file)
ds_b_b = NCDataset(weight_bool_basic_file)
ds_b_p = NCDataset(weight_bool_precise_file)

# Transformation des NCDatasets en matrice de dim 2, sur la var de temp
weights_prop_basic = ds_p_b["final_weights"][:,:]
# weights_prop_precise = ds_p_p["final_weights"][:,:]
weights_bool_basic = ds_b_b["mask"][:,:]
weights_bool_precise = ds_b_p["mask"][:,:]

# Fermeture des NCdataset pour éviter trop de poids sur la RAM
close(ds_p_b)
#close(ds_p_p)
close(ds_b_b)
close(ds_b_p)



"""
    compute_general_climatology(data_folder, weights, year_range; mode, ...)

Fonction principale (Orchestrateur) pour le calcul de climatologies.
Elle exécute séquentiellement :
1. L'accumulation des données brutes depuis les fichiers NetCDF.
2. Le calcul mathématique des moyennes et le masquage.
3. L'exportation optionnelle vers un nouveau fichier NetCDF.

Arguments :
- `data_folder` : Chemin vers le dossier des fichiers ERA5 bruts.
- `weights` : Matrice de poids spatiaux (utilisée pour créer le masque visuel).
- `year_range` : Plage d'années à traiter (ex: 1950:2020).
- `mode` : Niveau d'agrégation (:monthly, :yearly, ou :total).
- `selected_months` : Vecteur d'entiers pour filtrer les mois (ex: [6, 7, 8] pour l'été).
- `selected_days` : (Optionnel) Vecteur d'entiers pour filtrer des jours spécifiques.
- `export_path` : (Optionnel) Chemin de fichier pour sauvegarder le résultat NetCDF.

Retourne :
- `final_cube` : Une matrice 3D [Lon, Lat, Temps] contenant les températures moyennes en Celsius.
"""
function compute_general_climatology(
    data_folder::String, 
    weights::Matrix{Float64}, 
    year_range; 
    mode::Symbol=:total, 
    selected_months=collect(1:12), 
    selected_days=nothing, 
    selected_hours=0:23,
    variable_name="t2m",
    export_path=nothing
)    

    # Standardization of the month/year to vector for the loops
    if year_range isa Integer
        year_range = [year_range]
    end

    if selected_months isa Integer
        selected_months = [selected_months]
    end
    if selected_hours isa Interger
        selected_hours = [selected_hours]
    end
    if selected_days isa Integer
        selected_days = [selected_days]
    end
    # 1. Accumulate Raw Data
    println("Step 1: Accumulating data...")
    # Pass selected_hours to the function
    (sums, counts, lons, lats) = accumulate_data(
        data_folder, year_range, mode, selected_months, selected_days, selected_hours, variable_name)

    # 2. Compute Means & Cube 
    println("Step 2: Computing means...")
    (final_cube, valid_times) = finalize_cube(sums, counts, weights, mode)

    # 3. Export (Optional)
    if !isnothing(export_path)
        println("Step 3: Exporting...")
        save_climatology_netcdf(
            export_path, final_cube, valid_times, lons, lats, mode)
    end

    return final_cube
end


"""
    accumulate_data(data_folder, year_range, mode, months, days, var_name)

Étape 1 du pipeline : Lecture et Accumulation.
Cette fonction parcourt les fichiers NetCDF un par un et accumule les températures dans des dictionnaires, sans jamais charger toutes les données en mémoire vive.

Fonctionnement :
- Elle détermine une "Clé d'Agrégation" pour chaque image (ex: "Janvier 1950" ou juste "1950") selon le `mode`.
- Elle somme les températures valides dans `sums_dict`.
- Elle compte le nombre de pixels valides dans `counts_dict`.

Arguments :
- `data_folder` : Dossier source.
- `year_range` : Années à parcourir.
- `mode` : Définit la logique de création des clés (:monthly, :yearly, :total).
- `months`/`days` : Filtres temporels.
- `var_name` : Nom de la variable dans le NetCDF (par défaut "t2m").

Retourne un tuple :
- `(sums_dict, counts_dict, lons, lats)`
"""
function accumulate_data(data_folder, year_range, mode, months, days, hours, var_name)
    sums_dict = Dict{Any, Matrix{Float64}}()
    counts_dict = Dict{Any, Matrix{Int}}()
    lons, lats = nothing, nothing

    for year in year_range
        for month in months
            month_str = lpad(month, 2, '0')
            files = glob("*$(year)_$(month_str)*.nc", data_folder)
            
            # Check if file exists to avoid crash
            if isempty(files)
                continue
            end

            NCDataset(files[1]) do ds
                                # Capture coordinates once
                if isnothing(lons)
                    lons, lats = ds["longitude"][:], ds["latitude"][:]
                end
                # On réccupère le vecteur de date sur le mois
                times = ds["valid_time"][:]

                ### UPDATED: Robust filtering for Days AND Hours
                # We select index 'i' if:
                # 1. The day is in 'days' (or days is nothing)
                # 2. AND the hour is in 'hours'
                indices = findall(t -> 
                    (isnothing(days) || day(t) in days) && 
                    (hour(t) in hours), 
                    times
                )
                
                if !isempty(indices)
                    # Load only specific time steps
                    data_slice = ds[var_name][:, :, indices]
                    
                    for t in 1:length(indices)
                        time_val = times[indices[t]]
                        
                        # Get key based on new logic (including daily)
                        key = get_binning_key(mode, year, month, time_val)
                        # Test si la variable "days" est vide, si non on réccupère les indices dont le jour correspond.
                        if !haskey(sums_dict, key)
                            dims = size(data_slice)[1:2]
                            sums_dict[key] = zeros(Float64, dims)
                            counts_dict[key] = zeros(Int, dims)
                        end

                        frame = data_slice[:, :, t]
                        # Handling Missing/NaN
                        valid = .!ismissing.(frame) .& .!isnan.(frame)
                        if any(valid)
                            sums_dict[key][valid] .+= frame[valid]
                            counts_dict[key][valid] .+= 1
                        end
                    end
                end
            end
        end
        print("\rProcessed Year $year    ")
    end
    println("")
    return sums_dict, counts_dict, lons, lats
end




# Function d'aide pour l'utilisation des clées
function get_binning_key(mode, year, month, time_val)
    if mode == :monthly
        return Date(year, month, 1)
    elseif mode == :daily
        return Date(year, month, day(time_val))
    elseif mode == :yearly
        return year
    else
        return "Total"
    end
end


"""
    finalize_cube(sums_dict, counts_dict, weights, mode)

Étape 2 du pipeline : Finalisation Mathématique.
Transforme les données accumulées en un cube de données propre et calibré.

Opérations effectuées :
1. Calcul de la moyenne arithmétique : Moyenne = Somme / Compteur.
2. Application du Masque Visuel : Met à NaN les zones où `weights` < 0.9 (ex: Océans).
3. Conversion d'unités : Kelvin vers Celsius (T - 273.15).
4. Tri chronologique : Assure que les cartes sont dans l'ordre (Janvier, Février...).
5. Empilement : Transforme la liste de cartes 2D en une matrice 3D.

Arguments :
- `sums_dict` / `counts_dict` : Les dictionnaires remplis par `accumulate_data`.
- `weights` : La matrice de poids pour le masquage.
- `mode` : Nécessaire pour savoir comment trier les clés temporelles.

Retourne un tuple :
- `(final_3d_matrix, valid_times)`
"""
function finalize_cube(sums_dict, counts_dict, weights, mode)

    weights=weights'
    visual_mask = fill(NaN, size(weights))
    visual_mask[weights .> 0.0] .= 1.0

    # Sort keys chronologically
    all_keys = collect(keys(sums_dict))
    if mode != :total
        sort!(all_keys)
    end

    map_list = Matrix{Float64}[]
    valid_times = []

    for k in all_keys
        if any(counts_dict[k] .> 0)
            # Mean = Sum / Count
            mean_grid = sums_dict[k] ./ counts_dict[k]
            
            # Kelvin -> Celsius
            mean_grid .-= 273.15

            mean_grid .*= visual_mask

            
            push!(map_list, mean_grid)
            push!(valid_times, k)
        end
    end

    # Stack to 3D
    return cat(map_list..., dims=3), valid_times
end




"""
    save_climatology_netcdf(path, matrix, times, lons, lats, mode)

Étape 3 du pipeline : Exportation.
Crée un fichier NetCDF standardisé compatible avec les outils externes (QGIS, Python, Panoply).

Fonctionnalités :
- Définit dynamiquement la dimension temporelle ("time" pour mensuel, "year" pour annuel).
- Ajoute les métadonnées essentielles (_FillValue, units, long_name).
- Gère la suppression du fichier existant pour éviter les conflits d'écriture.

Arguments :
- `path` : Chemin de destination (ex: "output/climatology.nc").
- `matrix` : La matrice 3D calculée.
- `times` : Le vecteur d'axe temporel (Dates ou Entiers).
- `lons`/`lats` : Les vecteurs de coordonnées géographiques.
- `mode` : Influence le nommage des dimensions.
"""
function save_climatology_netcdf(path, matrix, times, lons, lats, mode)
    
    NCDataset(path, "c") do ds
        defDim(ds, "longitude", length(lons))
        defDim(ds, "latitude", length(lats))
        
        # Define time dimension name
        time_name = (mode == :yearly) ? "year" : "time"
        defDim(ds, time_name, length(times))

        defVar(ds, "longitude", Float64, ("longitude",))[:] = lons
        defVar(ds, "latitude", Float64, ("latitude",))[:] = lats
        
        # Write Time Variable
        if mode == :monthly || mode == :daily ### UPDATED
            # Both monthly and daily use Date objects
            defVar(ds, "time", times, ("time",))
        elseif mode == :yearly
            defVar(ds, "year", Int32, ("year",))[:] = times
        else
            defVar(ds, "time", [0.0], ("time",))
        end

        v_temp = defVar(ds, "temperature", Float64, ("longitude", "latitude", time_name), attrib = [
            "_FillValue" => NaN,
            "units" => "Celsius"])
        v_temp[:,:,:] = matrix
    end
end


matrice  = compute_general_climatology(data_folder_basic, weights_prop_basic, 2025; mode = :daily, selected_hours=18, selected_months=12)