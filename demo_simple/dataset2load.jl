# Assignation des différents chemins utiles
data_folder_precise = "data/france/france_9"
data_folder_basic ="data/france/france_31"
data_folder_ca_9="src/era5_land_ca_t2m"
data_folder_ca_31 = "src/ca_31km/era5_ca_t2m_31km"
weight_prop_basic = "data/masks/weights_france_31.nc"
weight_prop_precise = "data/masks/weights_france_9.nc"
weight_canada_9="data/masks/weights_canada_9.nc"
weight_canada_31="data/masks/weights_canada_31.nc"


# Création du dossier de sortie pour les plots
plot_dir = joinpath(@__DIR__, "plot")
isdir(plot_dir) || mkpath(plot_dir)