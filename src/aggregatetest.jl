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
- `global_file_path` : Chemin vers le fichier des données ERA5 bruts.
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
    global_file_path::String, 
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
    if selected_hours isa Integer
        selected_hours = [selected_hours]
    end
    if selected_days isa Integer
        selected_days = [selected_days]
    end
    # 1. Accumulate Raw Data
    println("Step 1: Accumulating data...")
    # Pass selected_hours to the function
    (sums, counts, lons, lats) = accumulate_data(
        global_file_path, year_range, mode, selected_months, selected_days, selected_hours, variable_name)

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
- `global_file_path` : Fichier source.
- `year_range` : Années à parcourir.
- `mode` : Définit la logique de création des clés (:monthly, :yearly, :total).
- `months`/`days` : Filtres temporels.
- `var_name` : Nom de la variable dans le NetCDF (par défaut "t2m").

Retourne un tuple :
- `(sums_dict, counts_dict, lons, lats)`
"""
function accumulate_data(global_file_path, year_range, mode, months, days, hours, var_name)
    # 1. Setup Thread-Safe Storage
    n_threads = Threads.nthreads()
    println("🚀 Starting Multithreaded Processing with $n_threads threads...")
    
    # Each thread gets its own dictionary to avoid conflicts (Race Conditions)
    thread_sums = [Dict{Any, Matrix{Float64}}() for _ in 1:n_threads]
    thread_counts = [Dict{Any, Matrix{Int}}() for _ in 1:n_threads]
    
    # 2. Get Metadata (Single Read)
    # We read dimensions once from the main thread
    lons, lats = nothing, nothing
    all_times = nothing
    
    NCDataset(global_file_path, "r") do ds
        time_name = haskey(ds, "valid_time") ? "valid_time" : "time"
        all_times = ds[time_name][:]
        lons = ds["longitude"][:]
        lats = ds["latitude"][:]
    end

    # 3. PARALLEL LOOP
    # Threads.@threads distributes the years among available cores 
    Threads.@threads for year in year_range
        tid = Threads.threadid()
        
        # Each thread opens its own file handle. 
        # This is safe for read-only and allows parallel decompression.
        NCDataset(global_file_path, "r") do ds
            
            # --- Fast Range Finder ---
            year_indices = findall(t -> Dates.year(t) == year, all_times)
            if isempty(year_indices) return end # Skip if empty

            range_start = minimum(year_indices)
            range_end = maximum(year_indices)
            
            # --- BULK READ (Decompression happens here) ---
            # Reading 1 year at once (~200MB) is efficient
            raw_year_data = ds[var_name][:, :, range_start:range_end]
            times_in_block = all_times[range_start:range_end]

            # --- Filter & Accumulate in RAM ---
            # We find indices relative to the loaded block (1 to 8760)
            relative_indices = findall(t -> 
                (Dates.month(t) in months) && 
                (isnothing(days) || Dates.day(t) in days) && 
                (isnothing(hours) || Dates.hour(t) in hours), 
                times_in_block
            )

            for i in relative_indices
                time_val = times_in_block[i]
                key = get_binning_key(mode, year, Dates.month(time_val), time_val)

                # Get the slice (View)
                frame = view(raw_year_data, :, :, i)
                
                # Check Dictionary for THIS thread
                if !haskey(thread_sums[tid], key)
                    dims = size(frame)
                    thread_sums[tid][key] = zeros(Float64, dims)
                    thread_counts[tid][key] = zeros(Int, dims)
                end

                # Vectorized Addition (Fast & Stable)
                valid_mask = .!ismissing.(frame) .& .!isnan.(frame)
                
                if any(valid_mask)
                    thread_sums[tid][key][valid_mask] .+= frame[valid_mask]
                    thread_counts[tid][key][valid_mask] .+= 1
                end
            end
        end
        # Print progress (Note: printing from threads can be messy, but useful)
        print(".") 
    end
    println("\n✅ Parallel processing complete. Merging results...")

    # 4. Merge results from all threads
    final_sums = merge_thread_dicts(thread_sums)
    final_counts = merge_thread_dicts(thread_counts)

    return final_sums, final_counts, lons, lats
end

# Helper to merge results from different threads
function merge_thread_dicts(dicts_list)
    merged = Dict{Any, Matrix{Float64}}()
    for d in dicts_list
        for (k, v) in d
            if haskey(merged, k)
                merged[k] .+= v
            else
                merged[k] = copy(v)
            end
        end
    end
    return merged
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


matrice  = compute_general_climatology("src/era5_global_1950_2025_basic.nc", weights_prop_basic, 1950:2025; mode = :yearly)