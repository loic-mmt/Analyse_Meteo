include("function2charge.jl")
include("dataset2load.jl")


# Exemple 1 : information sur une heure precise (03/10/1997 a 19h).
matrix1 = compute_general_climatology(data_folder_ca, weight_canada, 1997; selected_days=03, selected_months=10, selected_hours=19)
# p1 = carte instantanee de temperature pour cette date/heure.
vizumap(matrix1, weight_canada)
# Serie moyenne spatiale (ici un seul pas de temps car on a filtré une seule heure).
means_vector_calculation(matrix1, weight_prop_basic)


# Exemple 2 : carte moyenne 2020.
# Ici chaque pixel est moyenne sur toutes les dates/heures disponibles en 2020 (mode=:total).
matrix2 = compute_general_climatology(data_folder_ca, weight_canada, 1950, mode=:total)
# p2 = carte des moyennes temporelles par pixel en 2020.
vizumap(matrix2, weight_canada)
means_vector_calculation(matrix2, weight_canada)



# Exemple 3 : cartes horaires 2012.
matrix3 = compute_general_climatology(data_folder_ca, weight_canada, 2012, mode=:hourly)
vector3 = means_vector_calculation(matrix3, weight_canada)
# p3 = evolution temporelle de la moyenne spatiale, lissage sur 24 pas de temps.
trends_climate(vector3, window=24)
# Animation des cartes 2D au fil du temps (un frame par pas horaire).
animate_climatology(matrix3, weight_canada)

# Exemple 4 : tendances annuelles (1950->2025) sur la 2e quinzaine de janvier (jours 15 a 31).
matrix4 = compute_general_climatology(data_folder_basic, weight_prop_basic, 1950:2025; mode=:yearly, selected_months=01, selected_days=15:31)
vector4 = means_vector_calculation(matrix4, weight_prop_basic)
# p4 = courbe des moyennes annuelles sur la periode selectionnee (lissage fenetre 5).
p4 = trends_climate(vector4, window=5)
save_plot(p4, joinpath(plot_dir, "ex4_trends_window5.png"))
slope4, pvalue4 = calculate_trends_glm(matrix4, weight_prop_basic)
# p4_glm_005 = carte des pentes de tendance significatives au seuil p <= 0.05.
p4_glm_005 = glm_visu_trend(slope4, pvalue4)
save_plot(p4_glm_005, joinpath(plot_dir, "ex4_glm_p005.png"))
# p4_glm_015 = meme carte avec seuil p <= 0.15.
p4_glm_015 = glm_visu_trend(slope4, pvalue4, p_value=0.15)
save_plot(p4_glm_015, joinpath(plot_dir, "ex4_glm_p015.png"))


# Exemple 5 : tendance mensuelle (1950->2025).
matrix5 = compute_general_climatology(data_folder_basic, weight_prop_basic, 1950:2025, mode=:monthly)
vector5 = means_vector_calculation(matrix5, weight_prop_basic)
# p5 = evolution mensuelle de la moyenne spatiale, lissage fenetre 12.
p5 = trends_climate(vector5, window=12)
save_plot(p5, joinpath(plot_dir, "ex5_trends_window12.png"))

# Exemple 6 : tendance annuelle sur toute la periode (1950->2025), avec plusieurs variantes.
matrix6 = compute_general_climatology(data_folder_basic, weight_prop_basic, 1950:2025, mode=:yearly)
vector6 = means_vector_calculation(matrix6, weight_prop_basic)
# p6_w5 = courbe brute + moyenne glissante (fenetre 5).
p6_w5 = trends_climate(vector6, window=5)
save_plot(p6_w5, joinpath(plot_dir, "ex6_trends_window5.png"))
# p6_trend = courbe + droite de tendance lineaire globale.
p6_trend = trends_climate(vector6, trend=true)
save_plot(p6_trend, joinpath(plot_dir, "ex6_trends_trend.png"))
# p6_trend_cut = courbe + tendances separées avant/apres l'indice de coupure.
p6_trend_cut = trends_climate(vector6, trend=true, cutting=30)
save_plot(p6_trend_cut, joinpath(plot_dir, "ex6_trends_trend_cut30.png"))
