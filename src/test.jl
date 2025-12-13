using NCDatasets
using Glob
using Statistics
weight_file = "Analyse_Meteo/loic/weights_france_final.nc"
ds_w = NCDataset(weight_file)
year_ex = "Analyse_Meteo/data/raw-yearly-combined/era5_fr_t2m/era5_t2m_fr_1950_01.nc"
d0150 = NCDataset(year_ex)
time <- d0150["valid_time"][:,:]