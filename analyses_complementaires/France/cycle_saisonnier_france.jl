cd(joinpath(@__DIR__, "..", ".."))
include("../../demo_france/function2charge.jl")
include("../../demo_france/dataset2load.jl")

# Définition du dossier de sortie pour les graphiques
plot_dir = joinpath(@__DIR__, "plot")
isdir(plot_dir) || mkpath(plot_dir)


# PREPARATION DU MASQUE

function create_france_weights(data_folder, old_weights, new_weights)
    isfile(new_weights) && return new_weights
    
    println("Découpage du masque de poids aux dimensions exactes (62x43)...")
    
    data_file = glob("*.nc", data_folder)[1]
    ds_data = NCDataset(data_file)
    lon_data, lat_data = ds_data["longitude"][:], ds_data["latitude"][:]
    close(ds_data)
    
    ds_w = NCDataset(old_weights)
    lon_w, lat_w = ds_w["longitude"][:], ds_w["latitude"][:]
    weights_full = ds_w["final_weights"][:,:]
    close(ds_w)
    
    lon_idx = [argmin(abs.(lon_w .- lon)) for lon in lon_data]
    lat_idx = [argmin(abs.(lat_w .- lat)) for lat in lat_data]
    
    lon_idx = [argmin(abs.(lon_w .- lon)) for lon in lon_data]
    lat_idx = [argmin(abs.(lat_w .- lat)) for lat in lat_data]
    
    if size(weights_full) == (length(lon_w), length(lat_w))
        weights_cropped = weights_full[lon_idx, lat_idx]
    elseif size(weights_full) == (length(lat_w), length(lon_w))
        weights_cropped = permutedims(weights_full[lat_idx, lon_idx])
    else
        error("Dimensions inattendues pour le masque de poids : $(size(weights_full))")
    end
    
    NCDataset(new_weights, "c") do ds
        defDim(ds, "longitude", length(lon_data))
        defDim(ds, "latitude", length(lat_data))
        defVar(ds, "longitude", lon_data, ("longitude",))
        defVar(ds, "latitude", lat_data, ("latitude",))
        
        v1 = defVar(ds, "final_weights", Float64, ("longitude", "latitude"))
        v1[:,:] = weights_cropped
        v2 = defVar(ds, "weights_frac", Float64, ("longitude", "latitude"))
        v2[:,:] = weights_cropped
    end
    return new_weights
end

fichier_poids_france = joinpath(@__DIR__, "weights_france_exact.nc")
create_france_weights(data_folder_basic, weight_prop_basic, fichier_poids_france)


# REDEFINITION DES FONCTIONS

function finalize_cube(sums_dict, counts_dict, weights, mode)
    all_keys = collect(keys(sums_dict))
    if mode != :total
        sort!(all_keys)
    end

    valid_keys = filter(k -> any(counts_dict[k] .> 0), all_keys)
    
    if isempty(valid_keys)
        println("Warning: No valid data found.")
        return Array{Float64}(undef, 0, 0, 0), []
    end

    first_data = sums_dict[valid_keys[1]]
    (n_lon, n_lat) = size(first_data)
    n_time = length(valid_keys)

    if size(weights) == (n_lon, n_lat)
        aligned_weights = weights
        visual_mask = fill(NaN, (n_lon, n_lat))
        visual_mask[aligned_weights .> 0.0] .= 1.0
    elseif size(weights') == (n_lon, n_lat)
        aligned_weights = weights'
        visual_mask = fill(NaN, (n_lon, n_lat))
        visual_mask[aligned_weights .> 0.0] .= 1.0
    else
        visual_mask = ones(Float64, n_lon, n_lat) 
    end

    final_3d = zeros(Float64, n_lon, n_lat, n_time)
    valid_times = []

    for (i, k) in enumerate(valid_keys)
        mean_grid = sums_dict[k] ./ counts_dict[k]
        mean_grid .-= 273.15
        mean_grid .*= visual_mask
        final_3d[:, :, i] = mean_grid
        push!(valid_times, k)
    end

    return final_3d, valid_times
end


function means_vector_calculation(data_3d::AbstractArray{<:Union{Missing, Float64}, 3}, weights_file::String)
    temp_means = Float64[]

    ds=NCDataset(weights_file)
    if haskey(ds, "final_weights")
        weights = ds["final_weights"][:,:]
    else
        weights = ds["weights_frac"][:,:]
    end
    close(ds)
    
    (n_lon, n_lat, n_time) = size(data_3d)
    if size(weights) == (n_lon, n_lat)
        aligned_weights = weights
    elseif size(weights') == (n_lon, n_lat)
        aligned_weights = weights'
    else
        aligned_weights = ones(Float64, n_lon, n_lat) 
    end

    for t in 1:n_time
        map_slice = data_3d[:, :, t]
        weighted_map = map_slice .* aligned_weights
        valid_mask = .!ismissing.(weighted_map) .& .!isnan.(weighted_map)

        numerator = sum(weighted_map[valid_mask])
        denominator = sum(aligned_weights[valid_mask])
        
        mean_val = denominator > 0 ? numerator / denominator : NaN
        push!(temp_means, mean_val)
    end
    return temp_means
end


# ANALYSE ET GRAPHIQUES

# Comparaison des cycles (1950-1980 vs 1995-2025)
m_old = compute_general_climatology(data_folder_basic, fichier_poids_france, 1950:1980, mode=:monthly)
v_old = means_vector_calculation(m_old, fichier_poids_france)
cycle_old = mean(reshape(v_old, 12, :), dims=2)[:, 1]

m_recent = compute_general_climatology(data_folder_basic, fichier_poids_france, 1995:2025, mode=:monthly)
v_recent = means_vector_calculation(m_recent, fichier_poids_france)
cycle_recent = mean(reshape(v_recent, 12, :), dims=2)[:, 1]

p_comp = plot(1:12, cycle_old, label="1950-1980", linestyle=:dash, linewidth=3, color=:blue)
plot!(p_comp, 1:12, cycle_recent, label="1995-2025", linewidth=3, color=:red, 
      title="Décalage du cycle saisonnier en France", 
      xlabel="Mois",
      ylabel="Température Moyenne (°C)",
      xticks=(1:12, ["Jan", "Fév", "Mar", "Avr", "Mai", "Juin", "Juil", "Août", "Sep", "Oct", "Nov", "Déc"]),
      legend=:topleft,
      size=(800, 500))

save_plot(p_comp, joinpath(plot_dir, "comparaison_saisonniere_france.png"))


# Évolution des 4 saisons (1950-2025)
years_saison = 1950:2025

m_hiver = compute_general_climatology(data_folder_basic, fichier_poids_france, years_saison, mode=:yearly, selected_months=[12, 1, 2])
v_hiver = means_vector_calculation(m_hiver, fichier_poids_france)

m_printemps = compute_general_climatology(data_folder_basic, fichier_poids_france, years_saison, mode=:yearly, selected_months=[3, 4, 5])
v_printemps = means_vector_calculation(m_printemps, fichier_poids_france)

m_ete = compute_general_climatology(data_folder_basic, fichier_poids_france, years_saison, mode=:yearly, selected_months=[6, 7, 8])
v_ete = means_vector_calculation(m_ete, fichier_poids_france)

m_automne = compute_general_climatology(data_folder_basic, fichier_poids_france, years_saison, mode=:yearly, selected_months=[9, 10, 11])
v_automne = means_vector_calculation(m_automne, fichier_poids_france)

p1 = plot(years_saison, v_hiver, title="Hiver (DJF)", color=:blue, legend=false, linewidth=2)
p2 = plot(years_saison, v_printemps, title="Printemps (MAM)", color=:green, legend=false, linewidth=2)
p3 = plot(years_saison, v_ete, title="Été (JJA)", color=:red, legend=false, linewidth=2)
p4 = plot(years_saison, v_automne, title="Automne (SON)", color=:orange, legend=false, linewidth=2)

p_4saisons = plot(p1, p2, p3, p4, 
    layout=(2, 2), 
    size=(1000, 700),
    plot_title="Évolution des températures par saison en France (1950-2025)",
    xlabel="Années", 
    ylabel="Température (°C)",
    margin=5Plots.mm
)

save_plot(p_4saisons, joinpath(plot_dir, "evolution_4saisons_france.png"))