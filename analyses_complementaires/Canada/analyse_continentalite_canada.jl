# ANALYSE DE LA CONTINENTALITÉ (Points Côtiers vs Centraux)

# 1. Configuration des chemins et chargement
project_dir = abspath(joinpath(@__DIR__, "..", ".."))
include(joinpath(project_dir, "demo_canada", "function2charge.jl"))
include(joinpath(project_dir, "demo_canada", "dataset2load.jl"))

cd(project_dir)

# Dossier de sortie pour les graphiques
plot_dir = joinpath(@__DIR__, "plot")
isdir(plot_dir) || mkpath(plot_dir)


# 3. CALCUL DE LA CLIMATOLOGIE
println("Calcul de la climatologie mensuelle (1950-2025)...")
matrix_saison = compute_general_climatology(data_folder_ca_31, weight_canada_31, 1950:2025, mode=:monthly)

# 4. IDENTIFICATION DES POINTS 
mask_data = [all(isnan.(matrix_saison[i, j, :])) ? 0.0 : 1.0 for i in 1:size(matrix_saison,1), j in 1:size(matrix_saison,2)]

indices_terre = findall(mask_data .== 1.0)

if isempty(indices_terre)
    error("Aucune donnée de terre détectée dans matrix_saison.")
end

# Point Côtier : Le pixel de terre le plus à l'Ouest (Longitude min)
lons_terres = [idx[1] for idx in indices_terre]
idx_cotier = indices_terre[argmin(lons_terres)]

# Point Central : Le pixel de terre le plus proche du centre géographique
avg_lon = mean([idx[1] for idx in indices_terre])
avg_lat = mean([idx[2] for idx in indices_terre])
dist_centre = [abs(idx[1] - avg_lon) + abs(idx[2] - avg_lat) for idx in indices_terre]
idx_central = indices_terre[argmin(dist_centre)]

println("Points identifiés : Côtier $idx_cotier | Central $idx_central")

# 5. EXTRACTION ET CALCUL DU CYCLE MOYEN
# Extraction des séries temporelles
serie_cotier  = matrix_saison[idx_cotier[1], idx_cotier[2], :]
serie_central = matrix_saison[idx_central[1], idx_central[2], :]

# Repliement en 12 mois x N années et calcul de la moyenne par mois
cycle_cotier  = mean(reshape(serie_cotier, 12, :), dims=2)[:, 1]
cycle_central = mean(reshape(serie_central, 12, :), dims=2)[:, 1]

# Calcul des amplitudes pour les labels
amp_cotier  = round(maximum(cycle_cotier) - minimum(cycle_cotier), digits=1)
amp_central = round(maximum(cycle_central) - minimum(cycle_central), digits=1)

# 6. GRAPHIQUE FINAL
println("Étape 3 : Génération du graphique...")
p_cont = plot(1:12, cycle_cotier, 
    label="Côtier - Amplitude: $(amp_cotier)°C", 
    linewidth=4, color=:dodgerblue, marker=:circle, markersize=5)

plot!(p_cont, 1:12, cycle_central, 
    label="Central - Amplitude: $(amp_central)°C", 
    linewidth=4, color=:crimson, marker=:diamond, markersize=5,
    title="Analyse de Continentalité : Influence de l'Océan\n(Données ERA5 31km, Canada)",
    xlabel="Mois", ylabel="Température Moyenne (°C)",
    xticks=(1:12, ["Jan", "Fév", "Mar", "Avr", "Mai", "Juin", "Juil", "Août", "Sep", "Oct", "Nov", "Déc"]),
    legend=:bottomright, 
    grid=:xy,
    size=(900, 600)
)

save_plot(p_cont, joinpath(plot_dir, "comparaison_cotier_central_canada.png"))
