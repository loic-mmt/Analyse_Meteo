using NCDatasets # Pour le chargement des datasets .nc
using Glob #Pour pouvoir chercher des noms de fichiers dans un répèrtoire facillement.
using Statistics #Pour utiliser moyenne (mean) sur des matrices/dataset
using StatsPlots
using Dates # Pour utiliser le format "date" à l'intérieur des fichiers .nc
using GLM # Pour analyses de tendance (modèle linéaire simple)
using DataFrames 
using Base.Threads 
using RollingFunctions # Pour la création simple et rapide de vecteurs pondérés -> courbes pondérés
using ArchGDAL
using Proj



"""
    save_plot(plot_obj, plot_file)

Sauvegarde un objet graphique généré par `Plots.jl` dans un fichier.
Si un fichier du même nom existe déjà, il est supprimé avant l'écriture pour éviter les conflits.

# Arguments
- `plot_obj` : L'objet graphique à sauvegarder.
- `plot_file::String` : Le chemin et le nom du fichier de destination (ex: "output.png").
"""
function save_plot(plot_obj, plot_file)
    isfile(plot_file) && rm(plot_file)
    savefig(plot_obj, plot_file)
end




"""
    compute_general_climatology(data_folder, weights_file, year_range; kwargs...)

Fonction principale (Orchestrateur) pour le calcul de climatologies.
Elle exécute séquentiellement :
1. L'accumulation des données brutes depuis les fichiers NetCDF.
2. Le calcul mathématique des moyennes et le masquage.
3. L'exportation optionnelle vers un nouveau fichier NetCDF.

# Arguments
- `data_folder::String` : Chemin vers le dossier des fichiers ERA5 bruts.
- `weights_file::String` : Chemin vers le fichier NetCDF contenant la matrice de poids spatiaux.
- `year_range` : Plage d'années à traiter (ex: 1950:2020) ou entier unique.

# Mots-clés (kwargs)
- `mode::Symbol` : Niveau d'agrégation (:hourly, :daily, :monthly, :yearly, ou :total). Défaut: `:total`.
- `selected_months` : Vecteur d'entiers pour filtrer les mois (ex: [6, 7, 8] pour l'été). Défaut: `1:12`.
- `selected_days` : Vecteur d'entiers ou `nothing` pour filtrer des jours spécifiques.
- `selected_hours` : Plage horaire ou entier pour filtrer des heures spécifiques. Défaut: `0:23`.
- `variable_name::String` : Nom de la variable à extraire du NetCDF. Défaut: `"t2m"`.
- `export_path::String` : (Optionnel) Chemin pour sauvegarder le résultat NetCDF.

# Retourne
- `final_cube` : Une matrice 3D [Lon, Lat, Temps] contenant les températures moyennes en Celsius.
"""
function compute_general_climatology(
    data_folder::String, 
    weights_file::String, 
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
    ds=NCDataset(weights_file)
    weights=ds["final_weights"][:,:]
    close(ds)
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
    accumulate_data(data_folder, year_range, mode, months, days, hours, var_name)

Étape 1 du pipeline : Lecture et Accumulation.
Parcourt les fichiers NetCDF et accumule les températures valides en mémoire (streaming) sans surcharger la RAM.

# Logique
- Détermine une clé d'agrégation pour chaque pas de temps via `get_binning_key`.
- Ajoute les températures à `sums_dict` et incrémente le compteur de pixels dans `counts_dict`.

# Arguments
- `data_folder` : Dossier source.
- `year_range` : Années à parcourir.
- `mode` : Logique de création des clés temporelles.
- `months`/`days`/`hours` : Filtres temporels pour l'extraction.
- `var_name` : Nom de la variable géospatiale cible.

# Retourne
Un tuple `(sums_dict, counts_dict, lons, lats)`.
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



"""
    get_binning_key(mode, year, month, time_val)

Génère une clé d'agrégation dynamique pour grouper les données lors de l'accumulation.

# Arguments
- `mode::Symbol` : Le niveau de résolution temporelle désiré (:monthly, :daily, :hourly, :yearly).
- `year`, `month` : Entiers représentant l'année et le mois en cours de traitement.
- `time_val` : L'objet Date/DateTime exact extrait du fichier NetCDF.

# Retourne
Une clé (Date, entier, ou string) utilisée pour identifier la tranche temporelle dans le dictionnaire.
"""
function get_binning_key(mode, year, month, time_val)
    if mode == :monthly
        return Date(year, month, 1)
    elseif mode == :daily
        return Date(year, month, day(time_val))
    elseif mode == :hourly
        return time_val     
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

# Opérations
1. Calcul de la moyenne arithmétique (Somme / Compteur).
2. Conversion d'unités : Kelvin vers Celsius (T - 273.15).
3. Masquage Visuel : Applique `NaN` aux pixels où le poids est <= 0.0 (ex: océans).
4. Empilement : Construit une matrice 3D ordonnée chronologiquement.

# Arguments
- `sums_dict`, `counts_dict` : Dictionnaires générés par `accumulate_data`.
- `weights` : Matrice 2D des poids spatiaux.
- `mode` : Définit le type de tri chronologique à appliquer.

# Retourne
Un tuple `(final_3d_matrix, valid_times)`.
""" 
function finalize_cube(sums_dict, counts_dict, weights, mode)

    # 1. Prepare Mask
    weights = weights'
    visual_mask = fill(NaN, size(weights))
    visual_mask[weights .> 0.0] .= 1.0

    # 2. Sort keys chronologically
    all_keys = collect(keys(sums_dict))
    if mode != :total
        sort!(all_keys)
    end

    # 3. Filter valid keys first (to know the exact size)
    valid_keys = filter(k -> any(counts_dict[k] .> 0), all_keys)
    
    if isempty(valid_keys)
        println("⚠️ Warning: No valid data found.")
        return Array{Float64}(undef, 0, 0, 0), []
    end

    # 4. Pre-allocation
    # Get dimensions from the first valid entry
    first_data = sums_dict[valid_keys[1]]
    (n_lon, n_lat) = size(first_data)
    n_time = length(valid_keys)

    final_3d = zeros(Float64, n_lon, n_lat, n_time)
    valid_times = []

    # 5. Fill the Matrix loop (Safe & Fast)
    for (i, k) in enumerate(valid_keys)
        # Mean = Sum / Count
        mean_grid = sums_dict[k] ./ counts_dict[k]
        
        # Kelvin -> Celsius
        mean_grid .-= 273.15

        # Apply Mask
        mean_grid .*= visual_mask

        # Assign to the 3D Cube layer
        final_3d[:, :, i] = mean_grid
        
        push!(valid_times, k)
    end

    return final_3d, valid_times
end




"""
    save_climatology_netcdf(path, matrix, times, lons, lats, mode)

Étape 3 du pipeline : Exportation.
Crée un fichier NetCDF standardisé, incluant les dimensions et métadonnées nécessaires pour une utilisation dans des logiciels SIG (QGIS, Panoply).

# Arguments
- `path::String` : Chemin du fichier de destination.
- `matrix` : Matrice 3D de données (Longitude, Latitude, Temps).
- `times` : Vecteur représentant l'axe temporel.
- `lons`, `lats` : Vecteurs des coordonnées géographiques.
- `mode::Symbol` : Détermine le nommage de la dimension temporelle ("year" ou "time").
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
        if mode == :monthly || mode == :daily || mode == :hourly
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

"""
    means_vector_calculation(data_3d, weights_file)

Calcule une série temporelle de températures moyennes pondérées spatialement.
Réduit les dimensions spatiales (Longitude × Latitude) pour produire un vecteur 1D 
représentant l'évolution globale au cours du temps.

# Arguments
- `data_3d::AbstractArray{T, 3}` : Cube de données spatiales et temporelles.
- `weights_file::String` : Chemin vers le fichier NetCDF contenant la matrice des poids (frac).

# Retourne
- `Vector{Float64}` : Liste chronologique des moyennes pondérées calculées.
"""
function means_vector_calculation(
    data_3d::AbstractArray{<:Union{Missing, Float64}, 3}, 
    weights_file::String
)
    # Initialisation du vecteur
    temp_means = Float64[]

    # Pour avoir la même taille (latxlon)
    ds=NCDataset(weights_file)
    weights=ds["final_weights"][:,:]
    close(ds)
    weights = weights'


    # On utilise "size" pour trouver la taille de la variable temps (différents noms)
    (n_lon, n_lat, n_time) = size(data_3d)

    # Boucle sur la période de temps
    for t in 1:n_time
        # Découpage de la matrice pour ne prendre que la tranche de temps t (calculs plus rapides)
        map_slice = data_3d[:, :, t]
            
        # Application des poids
        weighted_map = map_slice .* weights
        
        # On garde seullement les non "missing" et "NaN"
        valid_mask = .!ismissing.(weighted_map) .& .!isnan.(weighted_map)

        # Sum(Data * Weight) / Sum(Weights)
        # On normalise pour avoir la moyenne de poids réelle
        numerator = sum(weighted_map[valid_mask])
        denominator = sum(weights[valid_mask])
        mean_val = numerator / denominator

        # On pousse dans le vecteur
        push!(temp_means, mean_val)
    end
    return temp_means
end

"""
    trends_climate(means; trend=false, cutting=Nothing, window=Nothing)

Trace l'évolution temporelle des températures et modélise les tendances.

# Arguments
- `means::Vector{Float64}` : Vecteur temporel des températures.

# Mots-clés
- `trend` : Si `true`, applique et trace une régression linéaire globale.
- `cutting` : Année (index) de rupture spatiale pour diviser la régression en deux périodes.
- `window` : Fenêtre temporelle (entier) pour superposer une courbe de moyenne mobile.

# Retourne
- `p` : L'objet Plot généré.
"""
function trends_climate(means::Vector{Float64}; trend = false, cutting=Nothing, window=Nothing)
    # Creation du dataframe pour les models et transfer en vecteur des années
    years_vec = Vector(1:length(means))
    df = DataFrame(Year = years_vec, Temp = means)

    # plot de base (scatter points)
    p = plot(df.Year, df.Temp,
        title = "Climate Trends Analysis",
        xlabel = "Year", ylabel = "Temperature (°C)",
        label = "Observed Means",
        seriestype = :scatter, 
        color = :pink,
        legend = :topleft,
        size = (800, 500),
        alpha=0.5
    )
    if window != Nothing
        roll_mean = runmean(means, window)
        plot!(p, years_vec, roll_mean,
            label = "$(window)-Moving Average", 
            color = :black, linewidth = 2)
    end

    
    if trend
        # Calcul du model et de la prédiction
        model_global = lm(@formula(Temp ~ Year), df)
        pred_global = predict(model_global, df)

        # ajout de la tendance sur le totalité
        slope_global = round(coef(model_global)[2], sigdigits= 2)
        confiance_global = round(coeftable(model_global).cols[6][2]-coeftable(model_global).cols[5][2], sigdigits=2)
        plot!(p, df.Year, pred_global,
            label = "Global trend, slope = $slope_global ± $confiance_global", color = :purple, linewidth=5)

        # boucle if si cutting est présent
        if cutting != Nothing
        
            # Tri des éléments plus petits que cutting, modélisation et plot
            df1 = filter(row -> row.Year <= cutting, df)
            model1 = lm(@formula(Temp ~ Year), df1)
            pred1 = predict(model1, df1)
            slope1 = round(coef(model1)[2], sigdigits= 2)
            confiance1 = round(coeftable(model1).cols[6][2]-coeftable(model1).cols[5][2], sigdigits=2)
            plot!(p, df1.Year, pred1, label = "Trend before $cutting , slope = $slope1 ± $confiance1", color = :yellow, linewidth=5)

            # Tri des éléments plus grands que cutting, modélisation et plot
            df2 = filter(row -> row.Year >= cutting, df)
            model2 = lm(@formula(Temp ~ Year), df2)
            pred2 = predict(model2, df2)
            slope2 = round(coef(model2)[2], sigdigits= 2)
            confiance2 = round(coeftable(model2).cols[6][2]-coeftable(model2).cols[5][2], sigdigits=2)
            plot!(p, df2.Year, pred2, label="Trend after $cutting , slope = $slope2 ± $confiance2", color = :orange, linewidth=5)
        end
    end

    display(p)
    return p
end

"""
    calculate_trends_glm(data_3d, weights_file)

Effectue une régression linéaire pixel par pixel sur le cube temporel complet.

Utilise un modèle GLM (`lm`) pour estimer la pente de température (°C/pas de temps) 
en chaque point géographique. Exclut les pixels masqués (poids < 0.5).

# Arguments
- `data_3d` : Cube de données [Lon, Lat, Temps].
- `weights_file::String` : Fichier contenant les poids pour ignorer les pixels invalides (océan).

# Retourne
- `slope_grid::Matrix{Float64}` : Carte 2D des coefficients directeurs (pentes).
- `p_value_grid::Matrix{Float64}` : Carte 2D des p-values associées à chaque pente.
"""
function calculate_trends_glm(data_3d::AbstractArray{<:Union{Missing, Float64}, 3}, weights_file)

    ds= NCDataset(weights_file)
    weights=ds["weights_frac"][:,:]
    close(ds)
    n_lon, n_lat, n_time = size(data_3d)
    
    # Initialize grids
    slope_grid = fill(NaN, n_lon, n_lat)
    p_value_grid = fill(NaN, n_lon, n_lat)
    
    time_steps = 1:n_time
    X_data = DataFrame(Time = time_steps)

    println("Calculating trends with GLM...")

    for i in 1:n_lon
        for j in 1:n_lat
            if weights[j, i] < 0.5
                continue
            end            
            y_temps = data_3d[i, j, :]
            

            if any(x -> ismissing(x) || isnan(x), y_temps)
                continue
            end
            X_data.Temp = y_temps
            
            model = lm(@formula(Temp ~ Time), X_data)
            
            # Extract Results
            # coef(model)[2] is the slope (Time coefficient)
            slope_grid[i, j] = coef(model)[2]
            
            # coeftable(model).cols[4][2] is the P-value for Time
            p_value_grid[i, j] = coeftable(model).cols[4][2]
        end
    end
    
    return slope_grid, p_value_grid
end

"""
    glm_visu_trend(slope_map, p_map; p_value=0.05, weights_file=nothing)

Génère une carte thermique (Heatmap) des tendances climatiques significatives.

Masque automatiquement les données non significatives où la p-value dépasse le seuil défini.

# Arguments
- `slope_map::AbstractArray` : Matrice 2D des pentes issues de `calculate_trends_glm`.
- `p_map::AbstractArray` : Matrice 2D des p-values.

# Mots-clés
- `p_value::Float64` : Seuil de significativité statistique (défaut 0.05).
- `weights_file` : (Optionnel) Fichier de poids utilisé pour projeter correctement les axes Longitude/Latitude.

# Retourne
- `p` : L'objet Plot de type Heatmap.
"""
function glm_visu_trend(slope_map::AbstractArray{Float64, 2}, p_map::AbstractArray{Float64, 2}; p_value=0.05, weights_file=nothing)

    # 2. Filter: Keep only significant trends (95% confidence)
    # We set non-significant pixels to NaN so they don't show up
    sig_slope_map = copy(slope_map)
    sig_slope_map[p_map .> p_value] .= NaN
    # 3. Plot
    limit = maximum(abs.(filter(!isnan, sig_slope_map)))
    lons = collect(1:size(sig_slope_map, 1))
    lats = collect(1:size(sig_slope_map, 2))
    if !isnothing(weights_file)
        ds = NCDataset(weights_file)
        lats = ds["latitude"][:]
        lons = ds["longitude"][:]
        close(ds)
    end
    z_map = sig_slope_map'
    idx_lon, idx_lat = sortperm(lons), sortperm(lats)
    x_plot, y_plot = lons[idx_lon], lats[idx_lat]
    z_plot = z_map[idx_lat, idx_lon]
    p = heatmap(
        x_plot,
        y_plot,
        z_plot,
        title = "Significant Warming Trends, p-value= $p_value",
        c = :balance,
        clims = (-limit, limit),
        #clims = (-0.065, 0.065),
        yflip = false,
        aspect_ratio = :equal
    )
    return p
end

"""
    animate_climatology(data_3d, weights_file; filename="temperature_evolution.gif")

Crée une animation GIF parcourant l'axe temporel du cube de données.
La limite des couleurs (`clims`) est verrouillée sur le minimum et maximum global 
de toute la période pour assurer la cohérence visuelle.

# Arguments
- `data_3d::AbstractArray` : Cube de températures [Lon, Lat, Temps].
- `weights_file::String` : Chemin vers le fichier NetCDF de masque/poids spatiaux.

# Mots-clés
- `filename::String` : Nom du fichier d'export de l'animation.

# Retourne
- Aucune valeur (Sauvegarde le fichier localement).
"""
function animate_climatology(data_3d::AbstractArray{<:Union{Missing, Float64}, 3}, weights_file; filename="temperature_evolution.gif")
    
    println("Generating animation...")

    # 1. Load coordinates and weights
    ds = NCDataset(weights_file)
    weights = ds["weights_frac"][:, :]
    lats = ds["latitude"][:]
    lons = ds["longitude"][:]
    close(ds)
    
    # --- Inline Axis Alignment Setup ---
    idx_lon, idx_lat = sortperm(lons), sortperm(lats)
    x_plot, y_plot = lons[idx_lon], lats[idx_lat]
    # -----------------------------------
    
    # 2. Determine fixed color limits for the whole period
    # We ignore NaNs so they don't break the min/max calculation
    valid_data = filter(!ismissing, data_3d)
    if isempty(valid_data)
        println("Error: Data contains only NaNs.")
        return
    end
    min_val, max_val = minimum(valid_data), maximum(valid_data)
    
    n_time = size(data_3d, 3)
    
    # 3. Create the Animation object
    anim = @animate for i in 1:n_time
        
        # Extract the 2D map for this year and transpose
        current_map = data_3d[:, :, i]'
        
        # --- Inline Masking ---
        current_map[weights .< 0.5] .= NaN
        
        # --- Apply axis alignment to the data ---
        z_plot = current_map[idx_lat, idx_lon]

        p = heatmap(
            x_plot,
            y_plot,
            z_plot,
            title = "Mean Temperature: $i",
            clims = (min_val, max_val),
            c = :thermal,   # Color palette
            xlabel = "Longitude",
            ylabel = "Latitude",
            aspect_ratio = :equal,
            right_margin = 5Plots.mm,
            yflip = false    # Give space for the colorbar
        )
        p
    end

    # 4. Save the GIF
    gif(anim, filename, fps = 5) 
    println("Saved animation to $filename")
end

"""
    vizumap(data_2d, weights_file)

Génère une carte thermique statique (Heatmap) à partir de données 2D ou de la première 
couche temporelle d'un cube 3D.

Applique un masque géographique via le fichier des poids.

# Arguments
- `data_2d` : Données spatiales en matrice (si 3D, sélectionne uniquement `[:, :, 1]`).
- `weights_file::String` : Fichier contenant la matrice des poids (frac) et les coordonnées lat/lon.

# Retourne
- `p` : L'objet Plot généré.
"""
function vizumap(data_2d, weights_file)
    ds = NCDataset(weights_file)
    weights = ds["weights_frac"][:,:]
    lats = ds["latitude"][:]
    lons = ds["longitude"][:]
    close(ds)
    
    valid_data = filter(!ismissing, data_2d)
    if isempty(valid_data)
        println("Error: Data contains only NaNs.")
        return
    end
    current_map = (ndims(data_2d) == 3) ? data_2d[:, :, 1]' : data_2d'
    
    # --- Inline Masking ---
    current_map[weights .< 0.5] .= NaN
    
    # --- Inline Axis Alignment ---
    idx_lon, idx_lat = sortperm(lons), sortperm(lats)
    x_plot, y_plot = lons[idx_lon], lats[idx_lat]
    z_plot = current_map[idx_lat, idx_lon]
    
    min_val, max_val = minimum(valid_data), maximum(valid_data)
    
    p = heatmap(
        x_plot,
        y_plot,
        z_plot,
        title = "Temperature",
        clims = (min_val, max_val),
        c = :thermal,
        xlabel = "Longitude",
        ylabel = "Latitude",
        aspect_ratio = :equal,
        right_margin = 5Plots.mm,
        yflip = false
    )
            
    return p
end
