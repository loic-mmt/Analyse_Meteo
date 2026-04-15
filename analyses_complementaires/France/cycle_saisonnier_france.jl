# ANALYSE SAISONNIÈRE ET ÉVOLUTION DES TEMPÉRATURES (FRANCE)

# 1. Configuration des chemins
project_dir = abspath(joinpath(@__DIR__, "..", ".."))
include(joinpath(project_dir, "demo_france", "function2charge.jl"))
include(joinpath(project_dir, "demo_france", "dataset2load.jl"))

cd(project_dir)

# Dossier de sortie pour les graphiques
plot_dir = joinpath(@__DIR__, "plot")
isdir(plot_dir) || mkpath(plot_dir)

# 2. REDÉFINITION DES FONCTIONS (Sécurité Dimensions et Kelvins)
function finalize_cube(sums_dict, counts_dict, weights, mode)
    all_keys = sort(collect(keys(sums_dict)))
    valid_keys = filter(k -> any(counts_dict[k] .> 0), all_keys)
    
    if isempty(valid_keys)
        return Array{Float64}(undef, 0, 0, 0), []
    end

    n_lon, n_lat = size(sums_dict[valid_keys[1]])
    
    # Alignement automatique du masque (transposition si nécessaire)
    if size(weights) == (n_lon, n_lat)
        aligned_weights = weights
    elseif size(weights') == (n_lon, n_lat)
        aligned_weights = weights'
    else
        aligned_weights = ones(Float64, n_lon, n_lat) 
    end

    visual_mask = [w > 0.0 ? 1.0 : NaN for w in aligned_weights]

    final_3d = zeros(Float64, n_lon, n_lat, length(valid_keys))
    for (i, k) in enumerate(valid_keys)
        mean_grid = (sums_dict[k] ./ counts_dict[k]) .- 273.15
        final_3d[:, :, i] = mean_grid .* visual_mask
    end
    return final_3d, valid_keys
end

function means_vector_calculation(data_3d::AbstractArray{<:Union{Missing, Float64}, 3}, weights_file::String)
    temp_means = Float64[]
    n_time = size(data_3d, 3)

    for t in 1:n_time
        # On extrait la tranche temporelle (la carte de France à l'instant t)
        map_slice = data_3d[:, :, t]
    
        valid_pixels = filter(!isnan, map_slice)
        
        if !isempty(valid_pixels)
            # On fait la moyenne arithmétique simple des pixels de terre
            push!(temp_means, mean(valid_pixels))
        else
            push!(temp_means, NaN)
        end
    end
    return temp_means
end

# 3. ANALYSE : COMPARAISON DES CYCLES (1950-1980 vs 1995-2025)
println("Calcul des cycles saisonniers...")

m_old = compute_general_climatology(data_folder_basic, weight_prop_basic, 1950:1985, mode=:monthly)
v_old = means_vector_calculation(m_old, weight_prop_basic)
cycle_old = mean(reshape(v_old, 12, :), dims=2)[:, 1]

m_recent = compute_general_climatology(data_folder_basic, weight_prop_basic, 1986:2025, mode=:monthly)
v_recent = means_vector_calculation(m_recent, weight_prop_basic)
cycle_recent = mean(reshape(v_recent, 12, :), dims=2)[:, 1]

p_comp = plot(1:12, cycle_old, label="1950-1985", ls=:dash, lw=3, color=:blue)
plot!(p_comp, 1:12, cycle_recent, label="1986-2025", lw=3, color=:red, 
      title="Réchauffement du cycle saisonnier en France", 
      ylabel="Température Moyenne (°C)",
      xticks=(1:12, ["Jan", "Fév", "Mar", "Avr", "Mai", "Juin", "Juil", "Août", "Sep", "Oct", "Nov", "Déc"]))

save_plot(p_comp, joinpath(plot_dir, "comparaison_saisonniere_france.png"))

# 4. ÉVOLUTION DES 4 SAISONS (1950-2025)
println("Calcul de l'évolution des saisons...")
years_saison = 1950:2025

# Extraction groupée par saison
saisons = [
    (name="Hiver", months=[12, 1, 2], col=:blue),
    (name="Printemps", months=[3, 4, 5], col=:green),
    (name="Été", months=[6, 7, 8], col=:red),
    (name="Automne", months=[9, 10, 11], col=:orange)
]

plots_saisons = []
for s in saisons
    m = compute_general_climatology(data_folder_basic, weight_prop_basic, years_saison, mode=:yearly, selected_months=s.months)
    v = means_vector_calculation(m, weight_prop_basic)
    p = plot(years_saison, v, title="$(s.name)", color=s.col, lw=2, legend=false)
    push!(plots_saisons, p)
end

p_4saisons = plot(plots_saisons..., layout=(2, 2), size=(1000, 700),
                  plot_title="Évolution des températures par saison (1950-2025)")

save_plot(p_4saisons, joinpath(plot_dir, "evolution_4saisons_france.png"))
