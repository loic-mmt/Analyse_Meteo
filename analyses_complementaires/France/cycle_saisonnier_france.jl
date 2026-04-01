cd(joinpath(@__DIR__, "..", ".."))
# Chargement des fonctions :
include("../../demo_france/function2charge.jl")
# Chargement des variables de base :
include("../../demo_france/dataset2load.jl")


# 1. Calculer la matrice mensuelle sur toute la période
matrix_saison = compute_general_climatology(data_folder_basic, weight_prop_basic, 1950:2025, mode=:monthly)
# 2. Calculer la moyenne spatiale (pour avoir une seule courbe pour toute la France)
vector_saison = means_vector_calculation(matrix_saison, weight_prop_basic)

# 3. Générer le graphique du cycle saisonnier
p_saison = plot(vector_saison, 
    title="Cycle saisonnier moyen en France (1950-2025)",
    xlabel="Mois", 
    ylabel="Température Moyenne (°C)",
    xticks=(1:12, ["Jan", "Fév", "Mar", "Avr", "Mai", "Juin", "Juil", "Août", "Sep", "Oct", "Nov", "Déc"]),
    linewidth=3,
    marker=:circle
)
save_plot(p_saison, joinpath(plot_dir, "cycle_saisonnier_total_france.png"))

# Comparer le cycle saisonnier avant et après

# Période ancienne (ex: 1950-1980)
m_old = compute_general_climatology(data_folder_basic, weight_prop_basic, 1950:1980, mode=:monthly)
v_old = means_vector_calculation(m_old, weight_prop_basic)

# Période récente (ex: 1995-2025)
m_recent = compute_general_climatology(data_folder_basic, weight_prop_basic, 1995:2025, mode=:monthly)
v_recent = means_vector_calculation(m_recent, weight_prop_basic)

# Graphique comparatif
p_comp = plot(v_old, label="1950-1980", linestyle=:dash, color=:blue)
plot!(v_recent, label="1995-2025", linewidth=2, color=:red, 
      title="Décalage du cycle saisonnier", 
      xticks=(1:12, ["J","F","M","A","M","J","J","A","S","O","N","D"]))

save_plot(p_comp, joinpath(plot_dir, "comparaison_saisonniere_france.png"))