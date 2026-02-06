# Chargement des fonctions :
include("function2charge.jl")
# Chargement des variables de base :
include("dataset2load.jl")



#Exemple 1 : Exctraction de l'information sur une heure précise : 
matrix1 = compute_general_climatology(data_folder_basic, weight_prop_basic, 1997; selected_days=03, selected_months=10, selected_hours=19)
vizumap(matrix1, weight_prop_basic)
means_vector_calculation(matrix1, weight_prop_basic)


# Exemple 2 : Extraction des informations de la carte des moyennes de température en France en 2020 :
matrix2 = compute_general_climatology(data_folder_basic, weight_prop_basic, 2020, mode=:total)
vizumap(matrix2, weight_prop_basic)
means_vector_calculation(matrix2, weight_prop_basic)

# Exemple 3 : Exctraction des cartes de toutes les heures de l'année 2012 : 
matrix3 = compute_general_climatology(data_folder_basic, weight_prop_basic, 2012, mode=:hourly)
vector3 = means_vector_calculation(matrix3, weight_prop_basic)
trends_climate(vector3, window=24)
animate_climatology(matrix3[:,:,24], weight_prop_basic)

#Exemple 4 : Extraction de tendance sur une période de temps (années 1950->2025) par années pour la deuxième quinzaine de janvier.
matrix4 = compute_general_climatology(data_folder_basic, weight_prop_basic, 1950:2025; mode=:yearly, selected_months=01, selected_days=15:31)
vector4 = means_vector_calculation(matrix4, weight_prop_basic)
trends_climate(vector4, window=5)
slope4, pvalue4 =calculate_trends_glm(matrix4, weight_prop_basic)
glm_visu_trend(slope4, pvalue4)
glm_visu_trend(slope4, pvalue4, p_value=0.15)


# Exemple 5 : Exctraction de tendance sur une période de temps (années 1980->2025) par mois.
matrix5 = compute_general_climatology(data_folder_basic, weight_prop_basic, 1950:2025, mode=:monthly)
vector5 = means_vector_calculation(matrix4, weight_prop_basic)
trends_climate(vector5, window=12)

#Exemple 6 : Exctraction de tendance sur la totalité de la période de temps (1950->2025) par ans.
matrix6 = compute_general_climatology(data_folder_basic, weight_prop_basic, 1950:2025, mode=:yearly)
vector6 = means_vector_calculation(matrix6, weight_prop_basic)
trends_climate(vector6, window=5)
trends_climate(vector6, trend=true)
trends_climate(vector6, trend=true, cutting=30)