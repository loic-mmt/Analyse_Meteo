# ANALYSE SAISONNIÈRE ET ÉVOLUTION DES TEMPÉRATURES (FRANCE)

# 1. Configuration des chemins
project_dir = abspath(joinpath(@__DIR__, "..", ".."))
include(joinpath(project_dir, "demo_france", "function2charge.jl"))
include(joinpath(project_dir, "demo_france", "dataset2load.jl"))

cd(project_dir)

# Dossier de sortie pour les graphiques
plot_dir = joinpath(@__DIR__, "plot")
isdir(plot_dir) || mkpath(plot_dir)


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


# 3. ÉVOLUTION DES ANOMALIES PAR SAISON (1950-2025)
println("Calcul des anomalies par saison...")
years_saison = 1950:2025
ref_period = 1960:1990
idx_ref = findall(y -> y in ref_period, years_saison)

saisons = [
    (name="Hiver", months=[12, 1, 2], col=:blue),
    (name="Printemps", months=[3, 4, 5], col=:green),
    (name="Été", months=[6, 7, 8], col=:red),
    (name="Automne", months=[9, 10, 11], col=:orange)
]

plots_saisons = []
for s in saisons
    m = compute_general_climatology(data_folder_basic, weight_prop_basic, years_saison, mode=:yearly, selected_months=s.months)
    v_abs = means_vector_calculation(m, weight_prop_basic)
    
    baseline = mean(v_abs[idx_ref])
    v_anom = v_abs .- baseline
    
    p = plot(years_saison, v_anom, 
             title="$(s.name)", color=s.col, lw=2, legend=false, ylabel="Anomalie (°C)")
    hline!(p, [0], color=:black, ls=:dash, lw=1.5)
    
    push!(plots_saisons, p)
end

p_4saisons = plot(plots_saisons..., layout=(2, 2), size=(1000, 700),
                  plot_title="Anomalies de température par saison (Réf: 1960-1990)")

save_plot(p_4saisons, joinpath(plot_dir, "courbe_anomalies_4saisons_france.png"))


# Courbe anomalie temperature
window = 5  # taille de lissage

plots_saisons = []
for s in saisons
    m = compute_general_climatology(data_folder_basic, weight_prop_basic,
                                   years_saison,
                                   mode=:yearly,
                                   selected_months=s.months)

    v = means_vector_calculation(m, weight_prop_basic)

    # 🔹 Indices période de référence
    idx_ref = findall(y -> y in 1960:1990, years_saison)

    # 🔹 Moyenne climatologique
    ref_mean = mean(v[idx_ref])

    # 🔹 Anomalies
    anomalies = v .- ref_mean

    # 🔹 Courbe de référence lissée (sur période 1960-1990 uniquement)
    ref_values = v[idx_ref]

    ref_smooth = [mean(ref_values[max(1,i-window):min(length(ref_values),i+window)])
                  for i in eachindex(ref_values)]

    # 🔹 Étendre la courbe sur toute la période
    ref_curve_full = fill(NaN, length(years_saison))
    ref_curve_full[idx_ref] .= ref_smooth

    # 🔹 Plot
    p = plot(years_saison, anomalies,
             color=s.col, lw=2, label="Anomalies")

    plot!(years_saison, ref_curve_full .- ref_mean,
          color=:black, lw=2, linestyle=:dash, label="Réf 1960-1990 (lissée)")

    push!(plots_saisons, p)
end

p_4saisons_2 = plot(plots_saisons..., layout=(2, 2), size=(1000, 700))

save_plot(p_4saisons_2, joinpath(plot_dir, "Anomalies_saisonniere_france.png"))
