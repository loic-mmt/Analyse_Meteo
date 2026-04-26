# Chargement des fonctions :
include("function2charge.jl")
# Chargement des variables de base :
include("dataset2load.jl")


# =====================================================================
# Exemple 1a : Information sur une heure precise (03/01/1997 a 19h) - FRANCE
# =====================================================================
println("\n--- Exemple 1a : France (03/01/1997 19:00) ---")
matrix1_fr = compute_general_climatology(data_folder_fr_31, weight_france_31, 1997; selected_days=3, selected_months=1, selected_hours=19)

p1_fr = vizumap(matrix1_fr, weight_france_31, "Température France - 3 janvier 1997 à 19h")
save_plot(p1_fr, joinpath(plot_dir, "ex1a_fr_vizumap.png"))

mean_1a = means_vector_calculation(matrix1_fr, weight_france_31)
println("Moyenne spatiale (France) : $(round(mean_1a[1], digits=2)) °C")

# =====================================================================
# Exemple 1b : Information sur une heure precise (03/01/1997 a 19h) - CANADA
# =====================================================================
println("\n--- Exemple 1b : Canada (03/01/1997 19:00) ---")
matrix1_ca = compute_general_climatology(data_folder_ca_31, weight_canada_31, 1997; selected_days=3, selected_months=1, selected_hours=19)

p1_ca = vizumap(matrix1_ca, weight_canada_31, "Température Canada - 3 janvier 1997 à 19h")
save_plot(p1_ca, joinpath(plot_dir, "ex1b_ca_vizumap.png"))

mean_1b = means_vector_calculation(matrix1_ca, weight_canada_31)
println("Moyenne spatiale (Canada) : $(round(mean_1b[1], digits=2)) °C")

# =====================================================================
# Exemple 2 : Carte moyenne annuelle 2020
# =====================================================================
println("\n--- Exemple 2 : Moyenne annuelle 2020 ---")
matrix2 = compute_general_climatology(data_folder_fr_31, weight_france_31, 2020, mode=:total)

p2 = vizumap(matrix2, weight_france_31, "Température moyenne annuelle France - 2020")
save_plot(p2, joinpath(plot_dir, "ex2_vizumap.png"))

mean_2 = means_vector_calculation(matrix2, weight_france_31)
println("Moyenne spatiale 2020 : $(round(mean_2[1], digits=2)) °C")

# ===== Canada =====
println("\n--- Exemple 2 Canada : Moyenne annuelle 2020 ---")
matrix2_ca = compute_general_climatology(data_folder_ca_31, weight_canada_31, 2020, mode=:total)
 
p2_ca = vizumap(matrix2_ca, weight_canada_31, "Température moyenne annuelle Canada - 2020")
save_plot(p2_ca, joinpath(plot_dir, "ex2_ca_vizumap.png"))
 
mean_2_ca = means_vector_calculation(matrix2_ca, weight_canada_31)
println("Moyenne spatiale Canada 2020 : $(round(mean_2_ca[1], digits=2)) °C")
 

# =====================================================================
# Exemple 3 : Cartes horaires 2012
# =====================================================================
println("\n--- Exemple 3 : Horaires 2012 ---")
matrix3 = compute_general_climatology(data_folder_fr_31, weight_france_31, 2012, mode=:hourly)
vector3 = means_vector_calculation(matrix3, weight_france_31)

p3 = trends_climate(vector3, window=24, title="Tendances horaires 2012 - France (fenêtre 24h)")
save_plot(p3, joinpath(plot_dir, "ex3_trends_window24.png"))


# ===== Canada =====
println("\n--- Exemple 3 : Horaires 2012 Canada ---")
matrix3_ca = compute_general_climatology(data_folder_ca_31, weight_canada_31, 2012, mode=:hourly)
vector3_ca = means_vector_calculation(matrix3_ca, weight_canada_31)
 
p3_ca = trends_climate(vector3_ca, window=24, title="Tendances horaires 2012 - Canada (fenêtre 24h)")
save_plot(p3_ca, joinpath(plot_dir, "ex3_ca_trends_window24.png"))
 

# =====================================================================
# Exemple 4 : Tendances annuelles (1950->2025) sur la 2e quinzaine de janvier
# =====================================================================
println("\n--- Exemple 4 : Tendances GLM (Janvier 15-31, 1950-2025) ---")
matrix4 = compute_general_climatology(data_folder_fr_31, weight_france_31, 1950:2025; mode=:yearly, selected_months=1, selected_days=15:31)
vector4 = means_vector_calculation(matrix4, weight_france_31)

p4 = trends_climate(vector4, window=5, title="Tendances annuelles France - Janvier (15-31) 1950-2025 (fenêtre 5 ans)")
save_plot(p4, joinpath(plot_dir, "ex4_trends_window5.png"))

slope4, pvalue4 = calculate_trends_glm(matrix4, weight_france_31)

p4_glm_005 = glm_visu_trend(slope4, pvalue4; weights_file=weight_france_31, title="Tendances GLM France - Janvier (15-31) p<0.05")
save_plot(p4_glm_005, joinpath(plot_dir, "ex4_glm_p005.png"))

p4_glm_015 = glm_visu_trend(slope4, pvalue4; p_value=0.15, weights_file=weight_france_31, title="Tendances GLM France - Janvier (15-31) p<0.15")
save_plot(p4_glm_015, joinpath(plot_dir, "ex4_glm_p015.png"))

matrix4 = nothing; slope4 = nothing; pvalue4 = nothing; GC.gc()

# ===== Canada =====
println("\n--- Exemple 4 : Tendances GLM (Janvier 15-31, 1950-2025) Canada ---")
matrix4_ca = compute_general_climatology(data_folder_ca_31, weight_canada_31, 1950:2025; mode=:yearly, selected_months=1, selected_days=15:31)
vector4_ca = means_vector_calculation(matrix4_ca, weight_canada_31)
 
p4_ca = trends_climate(vector4_ca, window=5, title="Tendances annuelles Canada - Janvier (15-31) 1950-2025 (fenêtre 5 ans)")
save_plot(p4_ca, joinpath(plot_dir, "ex4_ca_trends_window5.png"))
 
slope4_ca, pvalue4_ca = calculate_trends_glm(matrix4_ca, weight_canada_31)
 
p4_glm_005_ca = glm_visu_trend(slope4_ca, pvalue4_ca; weights_file=weight_canada_31, title="Tendances GLM Canada - Janvier (15-31) p<0.05")
save_plot(p4_glm_005_ca, joinpath(plot_dir, "ex4_ca_glm_p005.png"))
 
p4_glm_015_ca = glm_visu_trend(slope4_ca, pvalue4_ca; p_value=0.15, weights_file=weight_canada_31, title="Tendances GLM Canada - Janvier (15-31) p<0.15")
save_plot(p4_glm_015_ca, joinpath(plot_dir, "ex4_ca_glm_p015.png"))
matrix4_ca = nothing; slope4_ca = nothing; pvalue4_ca = nothing; GC.gc()


# =====================================================================
# Exemple 5 : Tendance mensuelle (1950->2025)
# =====================================================================
println("\n--- Exemple 5 : Evolution mensuelle (1950-2025) ---")
matrix5 = compute_general_climatology(data_folder_fr_31, weight_france_31, 1950:2025, mode=:monthly)
vector5 = means_vector_calculation(matrix5, weight_france_31)

p5 = trends_climate(vector5, window=12, title="Evolution mensuelle France 1950-2025 (fenêtre 12 mois)")
save_plot(p5, joinpath(plot_dir, "ex5_trends_window12.png"))

matrix5 = nothing; GC.gc()

# ===== Canada =====
println("\n--- Exemple 5 : Evolution mensuelle (1950-2025) Canada ---")
matrix5_ca = compute_general_climatology(data_folder_ca_31, weight_canada_31, 1950:2025, mode=:monthly)
vector5_ca = means_vector_calculation(matrix5_ca, weight_canada_31)
 
p5_ca = trends_climate(vector5_ca, window=12, title="Evolution mensuelle Canada 1950-2025 (fenêtre 12 mois)")
save_plot(p5_ca, joinpath(plot_dir, "ex5_ca_trends_window12.png"))
 
matrix5_ca = nothing; GC.gc()

# =====================================================================
# Exemple 6 : Tendance annuelle globale (1950->2025)
# =====================================================================
println("\n--- Exemple 6 : Regressions annuelles (1950-2025) ---")
matrix6 = compute_general_climatology(data_folder_fr_31, weight_france_31, 1950:2025, mode=:yearly)
vector6 = means_vector_calculation(matrix6, weight_france_31)

p6_w5 = trends_climate(vector6, trend=true, window=5, title="Tendances annuelles France 1950-2025 (fenêtre 5 ans)")
save_plot(p6_w5, joinpath(plot_dir, "ex6_trends_window5.png"))

p6_trend = trends_climate(vector6, trend=true, title="Tendances annuelles France 1950-2025 (régression linéaire)")
save_plot(p6_trend, joinpath(plot_dir, "ex6_trends_trend.png"))

p6_trend_cut = trends_climate(vector6, trend=true, cutting=30, title="Tendances annuelles France 1950-2025 (régression linéaire - cutoff 30 ans)")
save_plot(p6_trend_cut, joinpath(plot_dir, "ex6_trends_trend_cut30.png"))


# ===== Canada =====
println("\n--- Exemple 6 : Regressions annuelles (1950-2025) Canada ---")
matrix6_ca = compute_general_climatology(data_folder_ca_31, weight_canada_31, 1950:2025, mode=:yearly)
vector6_ca = means_vector_calculation(matrix6_ca, weight_canada_31)
 
p6_ca_w5 = trends_climate(vector6_ca, trend=true, window=5, title="Tendances annuelles Canada 1950-2025 (fenêtre 5 ans)")
save_plot(p6_ca_w5, joinpath(plot_dir, "ex6_ca_trends_window5.png"))
 
p6_ca_trend = trends_climate(vector6_ca, trend=true, title="Tendances annuelles Canada 1950-2025 (régression linéaire)")
save_plot(p6_ca_trend, joinpath(plot_dir, "ex6_ca_trends_trend.png"))
 
p6_ca_trend_cut = trends_climate(vector6_ca, trend=true, cutting=30, title="Tendances annuelles Canada 1950-2025 (régression linéaire - cutoff 30 ans)")
save_plot(p6_ca_trend_cut, joinpath(plot_dir, "ex6_ca_trends_trend_cut30.png"))
 
matrix6_ca = nothing; GC.gc()



# =====================================================================
# Exemple 7 : Comparaison Spatiale Haute vs Basse Resolution (France, Canicule Aout 2003)
# =====================================================================
println("\n--- Exemple 7 : Resolution 31km vs 9km (France, Aout 2003) ---")

# Calcul sur la grille standard (31km)
matrix7_31 = compute_general_climatology(data_folder_fr_31, weight_france_31, 2003; selected_months=8, selected_days=1:15, mode=:total)
p7_31 = vizumap(matrix7_31, weight_france_31, title="31km Température Moyenne France 1-15 Août 2003")
save_plot(p7_31, joinpath(plot_dir, "ex7_fr_2003_31km.png"))

# Calcul sur la grille haute resolution (9km)
matrix7_9 = compute_general_climatology(data_folder_fr_9, weight_france_9, 2003; selected_months=8, selected_days=1:15, mode=:total)
p7_9 = vizumap(matrix7_9, weight_france_9, title="9km Température Moyenne France 1-15 Août 2003")
save_plot(p7_9, joinpath(plot_dir, "ex7_fr_2003_9km.png"))

# Liberation memoire immediate
matrix7_31 = nothing; matrix7_9 = nothing; GC.gc()


# =====================================================================
# Exemple 8 : Analyse de Tendance GLM Haute Resolution (Canada 9km).
# =====================================================================
println("\n--- Exemple 8 : Tendances GLM Canada 9km 1950-2025) ---")

matrix8 = compute_general_climatology(data_folder_ca_9, weight_canada_9, 1950:2025; mode=:yearly)

# Calcul du modèle de régression linéaire sur la matrice haute résolution
slope8, pvalue8 = calculate_trends_glm(matrix8, weight_canada_9)

# Filtrage et visualisation des pentes
p8_glm = glm_visu_trend(slope8, pvalue8; weights_file=weight_canada_9, title="Tendances GLM canada 1950-2025")
save_plot(p8_glm, joinpath(plot_dir, "ex8_ca_glm_9km.png"))

matrix8 = nothing; slope8 = nothing; pvalue8 = nothing; GC.gc()

# Analyse de Tendance GLM Haute Resolution (France 9km).

matrix8b = compute_general_climatology(data_folder_fr_9, weight_france_9, 1950:2025; mode=:yearly)

# Calcul du modèle de régression linéaire sur la matrice haute résolution
slope8b, pvalue8b = calculate_trends_glm(matrix8b, weight_france_9)

# Filtrage et visualisation des pentes
p8b_glm = glm_visu_trend(slope8b, pvalue8b; weights_file=weight_france_9, title="Tendances GLM France 1950-2025")
save_plot(p8b_glm, joinpath(plot_dir, "ex8_fr_glm_9km_scale.png"))

matrix8b = nothing; slope8b = nothing; pvalue8b = nothing; GC.gc()

# =====================================================================
# Exemple 9 : Animation Haute Résolution Mensuelle (France 9km - 2022)
# =====================================================================
println("\n--- Exemple 9 : Animation Mensuelle France 9km (2022) ---")

# Extraction des moyennes mensuelles successives
matrix9 = compute_general_climatology(data_folder_fr_9, weight_france_9, 2022; mode=:monthly)

plot_file9 = joinpath(plot_dir, "ex9_fr_animation_9km.gif")
animate_climatology(matrix9, weight_france_9, filename=plot_file9, title = "Evolution Températures Mensuelles France 2022")

matrix9 = nothing; GC.gc()

# ===== Canada =====
println("\n--- Exemple 9 : Animation Mensuelle Canada 9km (2022) ---")
 
# Extraction des moyennes mensuelles successives
matrix9_ca = compute_general_climatology(data_folder_ca_9, weight_canada_9, 2022; mode=:monthly)
 
plot_file9_ca = joinpath(plot_dir, "ex9_ca_animation_9km.gif")
animate_climatology(matrix9_ca, weight_canada_9, filename=plot_file9_ca, title = "Evolution Températures Mensuelles Canada 2022")
 
matrix9_ca = nothing; GC.gc()

# =====================================================================
# Exemple 10 : Serie Temporelle Comparative (Canada été contre hiver)
# =====================================================================

println("\n--- Exemple 10 : Séries temporelles Canada (Hiver vs Ete, 1995-2025) ---")
 
# Extraction pour la période hivernale stricte (Dec, Jan, Fev)
matrix10_ca_winter = compute_general_climatology(data_folder_ca_31, weight_canada_31, 1995:2025; mode=:yearly, selected_months=[12, 1, 2])
vector10_ca_winter = means_vector_calculation(matrix10_ca_winter, weight_canada_31)
 
p10_ca_winter = trends_climate(vector10_ca_winter, trend=true, title="Tendances hivernales Canada 1995-2025")
save_plot(p10_ca_winter, joinpath(plot_dir, "ex10_ca_winter_trends_31km.png"))
 
# Extraction pour la période estivale stricte (Juin, Juil, Aout)
matrix10_ca_summer = compute_general_climatology(data_folder_ca_31, weight_canada_31, 1995:2025; mode=:yearly, selected_months=[6, 7, 8])
vector10_ca_summer = means_vector_calculation(matrix10_ca_summer, weight_canada_31)
 
p10_ca_summer = trends_climate(vector10_ca_summer, trend=true, title="Tendances estivales Canada 1995-2025")
save_plot(p10_ca_summer, joinpath(plot_dir, "ex10_ca_summer_trends_31km.png"))
 
matrix10_ca_winter = nothing; matrix10_ca_summer = nothing; GC.gc()

# =====================================================================
# Exemple 11 : Dynamique Horizontale Haute Résolution (Canada 31km)
# =====================================================================
println("\n--- Exemple 11 : Dynamique France 9km (Janvier 2018) ---")

# Génération d'une carte par heure pour observer la propagation spatiale des vagues de froids en janvier
matrix11 = compute_general_climatology(data_folder_fr_9, weight_france_9, 2018; mode=:hourly, selected_months=1, selected_hours=[0,6,12,18])

plot_file11 = joinpath(plot_dir, "ex11_fr_hourly_propagation_9km.gif")
animate_climatology(matrix11, weight_canada_31, filename=plot_file11, title="Evolution Températures Janvier 2018 France")

matrix11 = nothing; GC.gc()

# ===== Canada =====
println("\n--- Exemple 11 : Dynamique Canada 31km (Janvier 2018) ---")
 
# Génération d'une carte par heure pour observer la propagation spatiale des vagues de froids en janvier
matrix11_ca = compute_general_climatology(data_folder_ca_31, weight_canada_31, 2018; mode=:hourly, selected_months=1, selected_hours=[0,6,12,18])
 
plot_file11_ca = joinpath(plot_dir, "ex11_ca_hourly_propagation_31km.gif")
animate_climatology(matrix11_ca, weight_canada_31, filename=plot_file11_ca, title="Evolution Températures Janvier 2018 Canada")
 
matrix11_ca = nothing; GC.gc()
