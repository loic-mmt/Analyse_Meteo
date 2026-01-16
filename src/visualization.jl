using NCDatasets # Pour le chargement des datasets .nc
using Glob #Pour pouvoir chercher des noms de fichiers dans un répèrtoire facillement.
using Statistics #Pour utiliser moyenne (mean) sur des matrices/dataset
using StatsPlots
using Dates # Pour utiliser le format "date" à l'intérieur des fichiers .nc
using GLM # Pour analyses de tendance (modèle linéaire simple)
using DataFrames 
using Base.Threads 
using RollingFunctions # Pour la création simple et rapide de vecteurs pondérés -> courbes pondérés


MMB_file = "data/processed-means/mean_months_basic.nc"
MMP_file = "data/processed-means/mean_months_precise.nc"
MYB_file = "data/processed-means/mean_years_basic.nc"
MYP_file = "data/processed-means/mean_years_precise.nc"
MMB_nc = NCDataset(MMB_file)
MMP_nc = NCDataset(MMP_file)
MYB_nc = NCDataset(MYB_file)
MYP_nc = NCDataset(MYP_file)
MMB_matrix = MMB_nc["temperature"][:,:,:]
MMP_matrix = MMP_nc["temperature"][:,:,:]
MYB_matrix = MYB_nc["temperature"][:,:,:]
MYP_matrix = MYP_nc["temperature"][:,:,:]

weight_prop_basic_file = "data/masks/weights_prop_basic.nc"
weight_prop_precise_file = "data/masks/weights_prop_precise.nc"
weight_bool_basic_file = "data/masks/weights_bool_basic.nc"
weight_bool_precise_file = "data/masks/weights_bool_precise.nc"

# Importation des poids
ds_p_b = NCDataset(weight_prop_basic_file)
# ds_p_p = NCDataset(weight_prop_precise_file)
ds_b_b = NCDataset(weight_bool_basic_file)
ds_b_p = NCDataset(weight_bool_precise_file)

# Transformation des NCDatasets en matrice de dim 2, sur la var de temp
weights_prop_basic = ds_p_b["final_weights"][:,:]
# weights_prop_precise = ds_p_p["final_weights"][:,:]
weights_bool_basic = ds_b_b["mask"][:,:]
weights_bool_precise = ds_b_p["mask"][:,:]

# Fermeture des NCdataset pour éviter trop de poids sur la RAM
close(ds_p_b)
#close(ds_p_p)
close(ds_b_b)
close(ds_b_p)
close(MMB_nc)
close(MMP_nc)
close(MYB_nc)
close(MYP_nc)



function means_vector_calculation(
    data_3d::AbstractArray{<:Union{Missing, Float64}, 3}, 
    weights::Matrix{Float64}
)
    # Initialisation du vecteur
    temp_means = Float64[]

    # Pour avoir la même taille (latxlon)
    weights = weights'


    # On utilise "size" pour trouver la taille de la variable temps (différents noms)
    (n_lon, n_lat, n_time) = size(data_3d)

    # Boucle sur la période de temps
    for t in 1:n_time
        # Découpage de la matrice pour ne prendre que la tranche de temps t (calculs plus rapides)
        map_slice = data_3d[:, :, t]
            
        # Application des poids
        weighted_map = map_slice .* weights
        
        # On garde seullement les non "missing"
        valid_mask = .!ismissing.(weighted_map)

        # Sum(Data * Weight) / Sum(Weights)
        # On normalise pour avoir la moyenne de poids réelle
        numerator = sum(weighted_map[valid_mask])
        denominator = sum(weights[valid_mask])
        mean_val = numerator / denominator

        # On pousse dans le vecteur
        push!(temp_means, mean_val)
    end
    return temp_means
end

vect_months_prop = means_vector_calculation(MMB_matrix, weights_prop_basic)


