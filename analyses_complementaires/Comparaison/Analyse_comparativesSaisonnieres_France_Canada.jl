# Analyse comparative des anomalies saisonnières (France vs Canada)

# 1. Configuration des chemins et chargement des dépendances
project_dir = abspath(joinpath(@__DIR__, "..", ".."))

# Chargement des fonctions principales
include(joinpath(project_dir, "demo_france", "function2charge.jl"))

# Chargement des variables 
include(joinpath(project_dir, "demo_france", "dataset2load.jl"))

cd(project_dir)

# Création du dossier de sortie pour les graphiques si nécessaire
plot_dir = joinpath(@__DIR__, "plot")
isdir(plot_dir) || mkpath(plot_dir)

# 2. Configuration de l'analyse

println("Début du calcul des anomalies croisées France - Canada...")
years_saison = 1950:2025

# Définition de la période de référence (Normale climatique)
ref_period = 1960:1990
idx_ref = findall(y -> y in ref_period, years_saison)

# Définition des 4 saisons (mois correspondants)
saisons = [
    (name="Hiver", months=[12, 1, 2]),
    (name="Printemps", months=[3, 4, 5]),
    (name="Été", months=[6, 7, 8]),
    (name="Automne", months=[9, 10, 11])
]

plots_saisons_comp = []

# 3. Calcul et de création des graphiques par saison
for s in saisons
    println("Traitement de la saison : $(s.name)...")

    # --- CALCULS POUR LA FRANCE ---
    m_fr = compute_general_climatology(data_folder_basic, weight_prop_basic, years_saison, mode=:yearly, selected_months=s.months)
    v_abs_fr = means_vector_calculation(m_fr, weight_prop_basic)
    
    # Point zéro de la France
    baseline_fr = mean(v_abs_fr[idx_ref])
    v_anom_fr = v_abs_fr .- baseline_fr
    
    # --- CALCULS POUR LE CANADA ---
    m_ca = compute_general_climatology(data_folder_ca_31, weight_canada_31, years_saison, mode=:yearly, selected_months=s.months)
    v_abs_ca = means_vector_calculation(m_ca, weight_canada_31)
    
    # Point zéro du Canada
    baseline_ca = mean(v_abs_ca[idx_ref])
    v_anom_ca = v_abs_ca .- baseline_ca
    
    # --- CRÉATION DU GRAPHIQUE POUR LA SAISON ---
    p = plot(years_saison, v_anom_fr, 
             title="$(s.name)", 
             label="France",
             color=:blue, 
             lw=2, 
             ylabel="Anomalie (°C)",
             legend=:topleft)
             
    # 2. Ajout de la courbe Canada (ligne pointillée rouge)
    plot!(p, years_saison, v_anom_ca,
          label="Canada",
          color=:red,
          lw=2,
          ls=:dash)
             
    # 3. Ajout de la ligne de référence à 0 (Normale)
    hline!(p, [0], color=:black, ls=:dot, lw=1.5, label="")
    
    push!(plots_saisons_comp, p)
end


println("Génération du graphique final...")

# Assemblage des 4 graphiques dans une grille 2x2
p_comp_final = plot(plots_saisons_comp..., layout=(2, 2), size=(1200, 800),
                  plot_title="Comparaison des anomalies par saison : France vs Canada (Réf: 1960-1990)")

# Sauvegarde
output_file = joinpath(plot_dir, "anomalies_saisonnieres_comparees_FR_CA.png")
save_plot(p_comp_final, output_file)

println("Graphique sauvegardé avec succès dans : ", output_file)