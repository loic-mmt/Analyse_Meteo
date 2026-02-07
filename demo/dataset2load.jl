# Assignation des différents chemins utiles
data_folder_precise = "data/raw_monthly_combined/precise"
data_folder_basic ="data/raw_monthly_combined/basic"
file_weight_prop_basic = "data/masks/weights_prop_basic.nc"
file_weight_prop_precise = "data/masks/weights_prop_precise.nc"

# Importation des poids
ds_p_b = NCDataset(file_weight_prop_basic)
ds_p_p = NCDataset(file_weight_prop_precise)

# Transformation des NCDatasets en matrice de dim 2, sur la var de temp
weight_prop_basic = ds_p_b["final_weights"][:,:]
weight_prop_precise = ds_p_p["final_weights"][:,:]

# Fermeture des NCdataset pour éviter trop de poids sur la RAM
close(ds_p_b)
close(ds_p_p)

# Création du dossier de sortie pour les plots
plot_dir = joinpath(@__DIR__, "plot")
isdir(plot_dir) || mkpath(plot_dir)