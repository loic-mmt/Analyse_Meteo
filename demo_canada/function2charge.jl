using NCDatasets # Pour le chargement des datasets .nc
using Glob #Pour pouvoir chercher des noms de fichiers dans un répèrtoire facillement.
using Statistics #Pour utiliser moyenne (mean) sur des matrices/dataset
using StatsPlots
using Dates # Pour utiliser le format "date" à l'intérieur des fichiers .nc
using GLM # Pour analyses de tendance (modèle linéaire simple)
using DataFrames 
using Base.Threads 
using RollingFunctions # Pour la création simple et rapide de vecteurs pondérés -> courbes pondérés


function save_plot(plot_obj, plot_file)
    if plot_obj === nothing
        @warn "Skipping save: plot object is nothing for $plot_file"
        return false
    end
    isfile(plot_file) && rm(plot_file)
    savefig(plot_obj, plot_file)
    return true
end

function country_mask_for_map(weights::AbstractMatrix{<:Real}, target_size::Tuple{Int, Int})
    mask = weights .> 0.0
    if size(mask) == target_size
        return mask
    end
    mask_t = permutedims(mask)
    if size(mask_t) == target_size
        return mask_t
    end
    error("Incompatible sizes between weights $(size(weights)) and map $(target_size).")
end

function add_country_outline!(p, mask::AbstractMatrix{Bool}; linecolor=:black, linewidth=1.2)
    contour!(
        p,
        Float64.(mask);
        levels=[0.5],
        c=linecolor,
        linewidth=linewidth,
        label=false
    )
    return p
end




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
- `mode` : Niveau d'agrégation (:daily, :monthly, :yearly, ou :total).
- `selected_months` : (Optionnel) Vecteur d'entiers pour filtrer les mois (ex: [6, 7, 8] pour l'été).
- `selected_days` : (Optionnel) Vecteur d'entiers pour filtrer des jours spécifiques.
- `selected_hours` : Plage horaire à traiter pour filter des heures spécifiques.
- `export_path` : (Optionnel) Chemin de fichier pour sauvegarder le résultat NetCDF.

Retourne :
- `final_cube` : Une matrice 3D [Lon, Lat, Temps] contenant les températures moyennes en Celsius.

Exporte :
- un fichier .nc contenant la matrice  `final_cube`.
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
        data_folder, year_range, mode, selected_months, selected_days, selected_hours, variable_name)

    # 2. Compute Means & Cube 
    println("Step 2: Computing means...")
    (final_cube, valid_times) = finalize_cube(sums, counts, weights, mode)
    isempty(valid_times) && error(
        "No valid data found for filters: years=$(year_range), months=$(selected_months), " *
        "days=$(selected_days), hours=$(selected_hours)."
    )

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
            files = if isfile(data_folder)
                [data_folder]
            elseif isdir(data_folder)
                glob("*$(year)_$(month_str)*.nc", data_folder)
            else
                error("Invalid data source '$data_folder': expected a .nc file or a directory.")
            end
            
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
                time_var = if haskey(ds, "valid_time")
                    "valid_time"
                elseif haskey(ds, "time")
                    "time"
                else
                    error("No time coordinate found (expected 'valid_time' or 'time').")
                end
                times = ds[time_var][:]

                ### UPDATED: Robust filtering for Days AND Hours
                # We select index 'i' if:
                # 1. The day is in 'days' (or days is nothing)
                # 2. AND the hour is in 'hours'
                indices = findall(t -> 
                    (Dates.year(t) == year) &&
                    (Dates.month(t) == month) &&
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
                        key = get_binning_key(mode, Dates.year(time_val), Dates.month(time_val), time_val)
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
    means_vector_calculation(data_3d, weights)

Calcule une série temporelle de températures moyennes pondérées spatialement.

Cette fonction réduit les dimensions spatiales (Longitude × Latitude) pour produire un vecteur 1D
représentant l'évolution de la moyenne globale au cours du temps. Elle gère les données
en ajustant la somme des poids dynamiquement.

# Arguments
- `data_3d::AbstractArray{T, 3}` : Cube de données [Longitude, Latitude, Temps].
- `weights::Matrix{Float64}` : Matrice de poids 2D. Elle est transposée en interne pour correspondre aux dimensions.

# Retourne
- `Vector{Float64}` : La liste des températures moyennes pour chaque pas de temps.
"""
function means_vector_calculation(
    data_3d::AbstractArray{<:Union{Missing, Float64}, 3}, 
    weights::Matrix{Float64}
)
    # Initialisation du vecteur
    temp_means = Float64[]

    # Pour avoir la même taille (latxlon)
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
    trends_climate(means, years_range; trend=false, cutting=0)

Affiche l'évolution temporelle des températures et calcule les tendances linéaires.

Permet de visualiser les données brutes sous forme de scatter plot, et d'ajouter en option
une régression linéaire globale ou segmentée (avant/après une date de coupure).

# Arguments
- `means::Vector{Float64}` : Le vecteur des températures moyennes (sortie de `means_vector_calculation`).
- `trend::Bool` (keyword) : Si `true`, calcule et affiche les droites de régression.
- `cutting::Int` (keyword) : Année de rupture pour calculer deux tendances distinctes (ex: 1980).

# Retourne
- Un objet `Plot` contenant le graphique.
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
    calculate_trends_glm(data_3d, weights)

Effectue une régression linéaire pixel par pixel sur l'ensemble de la grille spatiale.

Utilise un modèle GLM (`lm`) pour estimer la pente de température (°C/an) en chaque point,
en excluant les zones masquées par les poids (ex: océans).

# Arguments
- `data_3d` : Cube de données [Lon, Lat, Time].
- `weights` : Matrice de poids utilisée comme masque (seuil à 0.5).

# Retourne
- `slope_grid::Matrix` : Carte des pentes (coefficients directeurs).
- `p_value_grid::Matrix` : Carte des p-values associées (significativité statistique).
"""
function calculate_trends_glm(data_3d::AbstractArray{<:Union{Missing, Float64}, 3}, weights)

    n_lon, n_lat, n_time = size(data_3d)
    
    # Initialize grids
    slope_grid = fill(NaN, n_lon, n_lat)
    p_value_grid = fill(NaN, n_lon, n_lat)
    
    time_steps = 1:n_time
    X_data = DataFrame(Time = time_steps)

    println("Calculating trends with GLM...")

    for i in 1:n_lon
        for j in 1:n_lat
            if weights[j, i] < 0.0
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
    glm_visu_trend(slope_map, p_map; p_value=0.05)

Affiche une carte thermique (heatmap) des tendances de température.

Filtre les pixels non significatifs en se basant sur la carte des p-values.
Les zones où p > `p_value` sont masquées (NaN) et n'apparaissent pas sur la carte.

# Arguments
- `slope_map` : Matrice des pentes calculées.
- `p_map` : Matrice des p-values correspondantes.
- `p_value` : Seuil de significativité (défaut 0.05 pour 95% de confiance).
"""
function glm_visu_trend(slope_map::AbstractArray{Float64, 2}, p_map::AbstractArray{Float64, 2}; p_value=0.05, weights=nothing)

    # 2. Filter: Keep only significant trends (95% confidence)
    # We set non-significant pixels to NaN so they don't show up
    sig_slope_map = copy(slope_map)
    print(length(sig_slope_map[p_map .> p_value]))
    sig_slope_map[p_map .> p_value] .= NaN
    # 3. Plot
    limit = maximum(abs.(filter(!isnan, sig_slope_map)))
    print(limit)
    p = heatmap(sig_slope_map', 
        title = "Significant Warming Trends, p-value= $p_value",
        c = :balance,
        #clims = (-limit, limit),
        clims = (-0.065, 0.065),
        yflip = true,
        aspect_ratio = :equal
    )
    if !isnothing(weights)
        mask = country_mask_for_map(weights, size(sig_slope_map'))
        add_country_outline!(p, mask)
    end
    return p
end

"""
    animate_climatology(data_3d, weights, valid_years; filename="...")

Génère une animation GIF montrant l'évolution des cartes de température année par année.

Applique un masque binaire basé sur les poids (terre/mer) pour ne visualiser que les
zones d'intérêt. L'échelle de couleur est fixée globalement pour permettre la comparaison.

# Arguments
- `data_3d` : Cube de données.
- `weights` : Matrice de masque (seuil à 0.5).
- `valid_years` : Vecteur des années correspondant à la dimension temporelle.
- `filename` : Nom du fichier de sortie (défaut "temperature_evolution.gif").
"""
function animate_climatology(data_3d::AbstractArray{<:Union{Missing, Float64}, 3}, weights; filename="temperature_evolution.gif")
    
    println("Generating animation...")

    # 1. Determine fixed color limits for the whole period
    # We ignore NaNs so they don't break the min/max calculation
    valid_data = filter(!ismissing, data_3d)
    if isempty(valid_data)
        println("Error: Data contains only NaNs.")
        return
    end
    min_val, max_val = minimum(valid_data), maximum(valid_data)
    mask = country_mask_for_map(weights, size(data_3d[:, :, 1]'))
    
    # 2. Create the Animation object
    anim = @animate for i in 1:length(data_3d[1,1,:])
        
        # Extract the 2D map for this year
        # Transpose (') is usually needed because Julia arrays are Col-Major
        # but heatmap expects [x, y]. 
        current_map = data_3d[:, :, i]'
        current_map[.!mask] .= NaN

        p = heatmap(current_map,
            title = "Mean Temperature: $i",
            clims = (min_val, max_val),
            c = :thermal,   # Color palette
            xlabel = "Longitude",
            ylabel = "Latitude",
            aspect_ratio = :equal,
            right_margin = 5Plots.mm,
            yflip = true    # Give space for the colorbar
        )
        add_country_outline!(p, mask)
        p
    end

    # 3. Save the GIF
    # fps = frames per second. 
    gif(anim, filename, fps = 5) 
    println("Saved animation to $filename")
end


function vizumap(data_2d, weights)
    valid_data = filter(!ismissing, data_2d)
    if isempty(valid_data)
        println("Error: Data contains only NaNs.")
        return
    end
    current_map = (ndims(data_2d) == 3) ? data_2d[:, :, 1]' : data_2d'
    mask = country_mask_for_map(weights, size(current_map))
    current_map[.!mask] .= NaN
    min_val, max_val = minimum(valid_data), maximum(valid_data)
    p = heatmap(current_map,
            title = "Temperature",
            clims = (min_val, max_val),
            c = :thermal,   # Color palette
            xlabel = "Longitude",
            ylabel = "Latitude",
            aspect_ratio = :equal,
            right_margin = 5Plots.mm,
            yflip = true)    # Give space for the colorbar
    add_country_outline!(p, mask)
    return p
end
