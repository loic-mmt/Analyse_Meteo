# ANALYSE SAISONNIÈRE ET ÉVOLUTION DES TEMPÉRATURES (CANADA)

# 1. Configuration des chemins
project_dir = abspath(joinpath(@__DIR__, "..", ".."))
include(joinpath(project_dir, "demo_canada", "function2charge.jl"))
include(joinpath(project_dir, "demo_canada", "dataset2load.jl"))

cd(project_dir)

# Dossier de sortie pour les graphiques
plot_dir = joinpath(@__DIR__, "plot")
isdir(plot_dir) || mkpath(plot_dir)


# 3. ANALYSE : COMPARAISON DES CYCLES (1950-1980 vs 1995-2025)
println("Calcul des cycles saisonniers...")

m_old = compute_general_climatology(data_folder_ca_31, weight_canada_31, 1950:1985, mode=:monthly)
v_old = means_vector_calculation(m_old, weight_canada_31)
cycle_old = mean(reshape(v_old, 12, :), dims=2)[:, 1]

m_recent = compute_general_climatology(data_folder_ca_31, weight_canada_31, 1986:2025, mode=:monthly)
v_recent = means_vector_calculation(m_recent, weight_canada_31)
cycle_recent = mean(reshape(v_recent, 12, :), dims=2)[:, 1]


p_comp = plot(1:12, cycle_old, label="1950-1985", ls=:dash, lw=3, color=:blue)
plot!(p_comp, 1:12, cycle_recent, label="1986-2025", lw=3, color=:red, 
      title="Réchauffement du cycle saisonnier au Canada", 
      ylabel="Température Moyenne (°C)",
      xticks=(1:12, ["Jan", "Fév", "Mar", "Avr", "Mai", "Juin", "Juil", "Août", "Sep", "Oct", "Nov", "Déc"]))

save_plot(p_comp, joinpath(plot_dir, "comparaison_saisonniere_canada.png"))


# Extraction groupée par saison
saisons = [
    (name="Hiver", months=[12, 1, 2], col=:blue),
    (name="Printemps", months=[3, 4, 5], col=:green),
    (name="Été", months=[6, 7, 8], col=:red),
    (name="Automne", months=[9, 10, 11], col=:orange)
]

plots_saisons = []
for s in saisons
    m = compute_general_climatology(data_folder_ca_31, weight_canada_31, 1950:2025, mode=:yearly, selected_months=s.months)
    v = means_vector_calculation(m, weight_canada_31)
    p = plot(1950:2025, v, title="$(s.name)", color=s.col, lw=2, legend=false)
    push!(plots_saisons, p)
end

p_4saisons = plot(plots_saisons..., layout=(2, 2), size=(1000, 700),
                  plot_title="Évolution des températures par saison (1950-2025)")

save_plot(p_4saisons, joinpath(plot_dir, "evolution_4saisons_canada.png"))
