# Assignation des différents chemins utiles
data_folder_fr_9 = "data/france/france_9"
data_folder_fr_31 ="data/france/france_31"
data_folder_ca_31="src/ca_31km"
data_folder_ca_9="src/ca_9km"
weight_france_31 = "data/masks/weights_france_31.nc"
weight_france_9 = "data/masks/weights_france_9.nc"
weight_canada_31="data/masks/weights_canada_31.nc"
weight_canada_9="data/masks/weights_canada_9.nc"

# Importation des poids
#ds_p_b = NCDataset(file_weight_prop_basic)
#ds_p_p = NCDataset(file_weight_prop_precise)
#ds_ca=NCDataset(file_ca_weights)

# Transformation des NCDatasets en matrice de dim 2, sur la var de temp
#weight_prop_basic = ds_p_b["final_weights"][:,:]
#weight_prop_precise = ds_p_p["final_weights"][:,:]
#weight_canada = ds_ca["final_weights"][:,:]

# Fermeture des NCdataset pour éviter trop de poids sur la RAM
#close(ds_p_b)
#close(ds_p_p)
#close(ds_ca)

# Création du dossier de sortie pour les plots
plot_dir = joinpath(@__DIR__, "plot")
isdir(plot_dir) || mkpath(plot_dir)