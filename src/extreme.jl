using NCDatasets
using Glob
using Statistics
using Plots
using Dates
using DataFrames
using CSV

# 1. CONFIGURATION ET CHARGEMENT DES DONNÉES

data_folder = "Analyse_Meteo/data/raw-yearly-combined/era5_fr_t2m"
weight_bool_file = "Analyse_Meteo/src/mask_france_boolean.nc"

# Chargement du Masque Binaire (France = 1, Mer/Ailleurs = 0)
println("Chargement du masque...")
ds_b = NCDataset(weight_bool_file)
mask_france = ds_b["mask"][:,:] # Matrice de 0 et 1
close(ds_b)

# Nombre total de pixels en France (pour faire la moyenne spatiale)
nb_pixels_france = sum(mask_france)

# 2. DÉFINITION DES SEUILS (Attention : ERA5 est en KELVIN)

# Seuil Canicule : 30°C
THRESHOLD_HOT_K = 30.0 + 273.15
# Seuil Gel : 0°C
THRESHOLD_COLD_K = 0.0 + 273.15

# 3. FONCTION DE CALCUL DES EXTRÊMES

function calculate_extremes(data_folder, mask, nb_pixels, year_range)

    # On va stocker les résultats ici
    results_hot = Float64[]
    results_cold = Float64[]
    years_list = Int[]

    println("Début de l'analyse des extrêmes sur $(year_range)...")

    for year in year_range
        # Compteurs pour l'année en cours
        total_hours_hot_france = 0
        total_hours_cold_france = 0

        # On cherche tous les fichiers de l'année (ex: 2024_01, 2024_02...)
        files = glob("*$(year)_*.nc", data_folder)

        if isempty(files)
            println("⚠️ Pas de données pour l'année $year")
            continue
        end

        # Pour chaque fichier (mois) de l'année
        for file in files
            NCDataset(file) do ds
                # Chargement de la variable t2m (Time, Lat, Lon)
                var = ds["t2m"]
                n_time = size(var)[3]

                # On boucle sur chaque HEURE du fichier
                for t in 1:n_time
                    # Lecture de la grille à l'instant t
                    temp_grid = var[:, :, t]

                    # CANICULE (> 30°C)

                    # On crée une carte de 1 (Chaud) et 0 (Pas chaud)
                    # On multiplie par le masque pour ne garder que la France
                    is_hot = (temp_grid .> THRESHOLD_HOT_K) .* mask

                    # On ajoute le nombre de pixels "chauds" à cet instant
                    total_hours_hot_france += sum(is_hot)

                    # ANALYSE GEL (< 0°C)

                    is_cold = (temp_grid .< THRESHOLD_COLD_K) .* mask
                    total_hours_cold_france += sum(is_cold)
                end
            end
        end

        # --- BILAN ANNUEL ---
        # On a la somme de toutes les heures chaudes sur tous les pixels.
        # Pour avoir "Le nombre d'heures moyen par pixel en France", on divise par le nb de pixels.
        avg_hours_hot = total_hours_hot_france / nb_pixels
        avg_hours_cold = total_hours_cold_france / nb_pixels

        push!(results_hot, avg_hours_hot)
        push!(results_cold, avg_hours_cold)
        push!(years_list, year)

        println("Année $year : $avg_hours_hot heures > 30°C (moyenne France)")
    end

    return years_list, results_hot, results_cold
end

# 4. EXÉCUTION ET SAUVEGARDE

years_to_process = 1950:1990

years, hot_hours, cold_hours = calculate_extremes(data_folder, mask_france, nb_pixels_france, years_to_process)

# --- VISUALISATION ---

# Graphique 1 : Canicules
p1 = plot(years, hot_hours,
    title="Évolution des heures de chaleur (> 30°C)",
    label="Heures > 30°C",
    color=:red,
    lw=2,
    xlabel="Année",
    ylabel="Nombre d'heures (Moyenne France)",
    legend=:topleft
)
savefig(p1, "evolution_canicule.png")
display(p1)

# Graphique 2 : Gel
p2 = plot(years, cold_hours,
    title="Évolution des heures de gel (< 0°C)",
    label="Heures < 0°C",
    color=:blue,
    lw=2,
    xlabel="Année",
    ylabel="Nombre d'heures (Moyenne France)",
    legend=:topright
)
savefig(p2, "evolution_gel.png")
display(p2)
