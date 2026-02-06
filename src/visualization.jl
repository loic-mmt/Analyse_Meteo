using NCDatasets # Pour le chargement des datasets .nc
using Glob #Pour pouvoir chercher des noms de fichiers dans un répèrtoire facillement.
using Statistics #Pour utiliser moyenne (mean) sur des matrices/dataset
using StatsPlots
using Dates # Pour utiliser le format "date" à l'intérieur des fichiers .nc
using GLM # Pour analyses de tendance (modèle linéaire simple)
using DataFrames 
using Base.Threads 
using RollingFunctions # Pour la création simple et rapide de vecteurs pondérés -> courbes pondérés


MMB_file = "data/processed-means/mean_months_basic.nc"
MMP_file = "data/processed-means/mean_months_precise.nc"
MYB_file = "data/processed-means/mean_years_basic.nc"
MYP_file = "data/processed-means/mean_years_precise.nc"
MMB_nc = NCDataset(MMB_file)
MMP_nc = NCDataset(MMP_file)
MYB_nc = NCDataset(MYB_file)
MYP_nc = NCDataset(MYP_file)
MMB_matrix = MMB_nc["temperature"][:,:,:]
MMP_matrix = MMP_nc["temperature"][:,:,:]
MYB_matrix = MYB_nc["temperature"][:,:,:]
MYP_matrix = MYP_nc["temperature"][:,:,:]

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


test_file = "output/climatology_yearly.nc"
testnc = NCDataset(test_file)
ttt = testnc["t2m"][:,:,:]
close(testnc)

ari =NCDataset("output/climatology_yearly_original.nc")
arfinal = ari["temperature"][:,:,:]


# Fermeture des NCdataset pour éviter trop de poids sur la RAM
close(ds_p_b)
#close(ds_p_p)
close(ds_b_b)
close(ds_b_p)
close(MMB_nc)
close(MMP_nc)
close(MYB_nc)
close(MYP_nc)


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
        
        # On garde seullement les non "missing"
        valid_mask = .!ismissing.(weighted_map)

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

vect_months_basic = means_vector_calculation(MMB_matrix, weights_prop_basic)
vect_years_basic = means_vector_calculation(MYB_matrix, weights_prop_basic)
vect_months_precise = means_vector_calculation(MM_matrix, weights_prop_precise)
vect_years_precise = means_vector_calculation(MYB_matrix, weights_prop_precise)

vect_test = means_vector_calculation(test_matrix, weights_prop_basic)
d_test = NCDataset("monthtest4.nc")
test_matrix = d_test["temperature"][:,:,:]

CDOtest = means_vector_calculation(ttt, weights_prop_basic)
aritest = means_vector_calculation(arfinal, weights_prop_basic)

trends_climate(vect_years_basic, 1950:2025)

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
function trends_climate(means::Vector{Float64}; trend = false, cutting=Nothing)
    # Creation du dataframe pour les models et transfer en vecteur des années
    years_vec = Vector(1:lenght(means))
    df = DataFrame(Year = years_vec, Temp = means)

    # plot de base (scatter points)
    p = plot(df.Year, df.Temp,
        title = "Climate Trends Analysis",
        xlabel = "Year", ylabel = "Temperature (°C)",
        label = "Observed Means",
        seriestype = :scatter, 
        color = :pink, alpha = 0.5,
        legend = :topleft,
        size = (800, 500)
    )
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
trends_climate(vect_test, 1950:2025, trend = true)
trends_climate(CDOtest, 1950:2025, trend = true)



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
    glm_visu_trend(slope_map, p_map; p_value=0.05)

Affiche une carte thermique (heatmap) des tendances de température.

Filtre les pixels non significatifs en se basant sur la carte des p-values.
Les zones où p > `p_value` sont masquées (NaN) et n'apparaissent pas sur la carte.

# Arguments
- `slope_map` : Matrice des pentes calculées.
- `p_map` : Matrice des p-values correspondantes.
- `p_value` : Seuil de significativité (défaut 0.05 pour 95% de confiance).
"""
function glm_visu_trend(slope_map::AbstractArray{Float64, 2}, p_map::AbstractArray{Float64, 2}; p_value=0.05)

    # 2. Filter: Keep only significant trends (95% confidence)
    # We set non-significant pixels to NaN so they don't show up
    sig_slope_map = copy(slope_map)
    print(length(sig_slope_map[p_map .> p_value]))
    sig_slope_map[p_map .> p_value] .= NaN
    # 3. Plot
    limit = maximum(abs.(filter(!isnan, sig_slope_map)))
    print(limit)
    heatmap(sig_slope_map', 
        title = "Significant Warming Trends",
        c = :balance,
        #clims = (-limit, limit),
        clims = (-0.065, 0.065),
        yflip = true,
        aspect_ratio = :equal
    )
end

glm_visu_trend(slope, p_value)


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
function animate_climatology(data_3d::AbstractArray{<:Union{Missing, Float64}, 3}, weights, valid_years::AbstractVector; filename="temperature_evolution.gif")
    
    println("Generating animation...")

    # 1. Determine fixed color limits for the whole period
    # We ignore NaNs so they don't break the min/max calculation
    valid_data = filter(!ismissing, data_3d)
    if isempty(valid_data)
        println("Error: Data contains only NaNs.")
        return
    end
    min_val, max_val = minimum(valid_data), maximum(valid_data)
    
    # 2. Create the Animation object
    anim = @animate for i in 1:length(valid_years)
        year = valid_years[i]
        
        # Extract the 2D map for this year
        # Transpose (') is usually needed because Julia arrays are Col-Major
        # but heatmap expects [x, y]. 
        current_map = data_3d[:, :, i]'
        current_map[weights .<0.5] .= NaN

        heatmap(current_map,
            title = "Mean Temperature: $year",
            clims = (min_val, max_val),
            c = :thermal,   # Color palette
            xlabel = "Longitude",
            ylabel = "Latitude",
            aspect_ratio = :equal,
            right_margin = 5Plots.mm,
            yflip = true    # Give space for the colorbar
        )
    end

    # 3. Save the GIF
    # fps = frames per second. 
    gif(anim, filename, fps = 5) 
    println("Saved animation to $filename")
end


"""
    animate_climatology(data_3d, weights, valid_times; filename="...")

Génère une animation GIF montrant l'évolution des cartes de température pas de temps par pas de temps.

Applique un masque binaire basé sur les poids (terre/mer) pour ne visualiser que les
zones d'intérêt.

# Arguments
- `data_3d` : Cube de données [Lon, Lat, Time].
- `weights` : Matrice de masque (seuil à 0.5).
- `valid_times` : Vecteur correspondant à la dimension temporelle (Années ou Dates).
- `filename` : Nom du fichier de sortie (défaut "temperature_evolution.gif").
"""
function animate_climatology(data_3d::AbstractArray{<:Union{Missing, Float64}, 3}, weights, valid_times::AbstractVector; filename="temperature_evolution.gif")
    
    println("Generating animation...")

    # 1. Determine fixed color limits for the whole period to allow comparison
    valid_data = filter(!ismissing, data_3d)
    if isempty(valid_data)
        println("Error: Data contains only NaNs.")
        return
    end
    min_val, max_val = minimum(valid_data), maximum(valid_data)
    
    # 2. Create the Animation object
    anim = @animate for i in 1:length(valid_times)
        # Fix: Use 'current_time' instead of 'year' to be generic (works for Date or Int)
        current_time = valid_times[i]
        
        # Extract the 2D map for this time step
        # Transpose (') to switch from [Lon, Lat] (Julia) to [x, y] (Plots.jl)
        current_map = data_3d[:, :, i]'
        
        # Apply visual mask
        current_map[weights .< 0.5] .= NaN

        heatmap(current_map,
            title = "Mean Temperature: $current_time", # Displays Year (Int) or Date (yyyy-mm-dd) automatically
            clims = (min_val, max_val),
            c = :thermal,   # Color palette
            xlabel = "Longitude",
            ylabel = "Latitude",
            aspect_ratio = :equal,
            right_margin = 5Plots.mm,
            yflip = true    # Orient North correctly
        )
    end

    # 3. Save the GIF
    gif(anim, filename, fps = 5) 
    println("✅ Saved animation to $filename")
end