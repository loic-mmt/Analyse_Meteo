# 1. Configuration des chemins et chargement
cd(joinpath(@__DIR__, "..", ".."))
include("../../demo_france/function2charge.jl")
include("../../demo_france/dataset2load.jl")

# Dossier de sortie spécifique
plot_dir = joinpath(@__DIR__, "plot")
isdir(plot_dir) || mkpath(plot_dir)

# Chargement du mask 
fichier_poids_france = joinpath(@__DIR__, "weights_france_exact.nc")

# Matrice mensuelle globale
matrix_saison = compute_general_climatology(data_folder_basic, fichier_poids_france, 1950:2025, mode=:monthly)

# Identification des points (Côtier vs Central)
ds_w = NCDataset(fichier_poids_france)
weights_map = ds_w["final_weights"][:,:]
close(ds_w)

# On récupère tous les indices de terre 
indices_terre = findall(weights_map .> 0.8)

# Point Côtier (Le plus à l'Ouest du masque de terre -> Bretagne)
idx_cotier = indices_terre[argmin([idx[1] for idx in indices_terre])]

# Point Central (calcul pour trouver le pixel qui se situe pile au centre géométrique de la France)
# On prend la moyenne des longitudes et latitudes valides
mid_lon = round(Int, mean([idx[1] for idx in indices_terre]))
mid_lat = round(Int, mean([idx[2] for idx in indices_terre]))
idx_central = (mid_lon, mid_lat)

# 3. Extraction et calcul du cycle moyen : idx_cotier[1] → position en longitude
# et idx_cotier[2] → position en latitude
serie_cotier  = matrix_saison[idx_cotier[1], idx_cotier[2], :]
serie_central = matrix_saison[idx_central[1], idx_central[2], :]

cycle_cotier  = mean(reshape(serie_cotier, 12, :), dims=2)[:, 1] 
# On transforme la série en tableau avec 12 lignes = les mois (janvier → décembre)
# et en colonnes = les années, puis on calcule la moyenne pour chaque mois
cycle_central = mean(reshape(serie_central, 12, :), dims=2)[:, 1]

# 4. Calcul de l'Amplitude Thermique (Indice de continentalité)
amp_cotier  = round(maximum(cycle_cotier) - minimum(cycle_cotier), digits=1)
amp_central = round(maximum(cycle_central) - minimum(cycle_central), digits=1)

# 5. Graphique
p_cont = plot(1:12, cycle_cotier, 
    label="Point Côtier (Bretagne) - Amp: $amp_cotier°C", 
    linewidth=4, color=:dodgerblue, marker=:circle)

plot!(p_cont, 1:12, cycle_central, 
    label="Point Central (Berry) - Amp: $amp_central°C", 
    linewidth=4, color=:orange, marker=:diamond,
    title="Continentalité : Impact de la proximité océanique",
    xlabel="Mois", ylabel="Température Moyenne (°C)",
    xticks=(1:12, ["Jan", "Fév", "Mar", "Avr", "Mai", "Juin", "Juil", "Août", "Sep", "Oct", "Nov", "Déc"]),
    legend=:bottomright,
    size=(900, 600)
)

save_plot(p_cont, joinpath(plot_dir, "comparaison_cotier_central.png"))