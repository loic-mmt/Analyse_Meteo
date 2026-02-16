include("function2charge.jl")
include("dataset2load.jl")


test_ca = NCDataset("demo/output/era5_ca_t2m/era5_t2m_ca_1950_02.nc")
matrix_ca = test_ca["t2m"][:,:,:]

test_original =NCDataset("data/canada/data_0.nc")
test_matrix=test_original["t2m"][:,:,:]

# Exemple 2 : carte moyenne 2020.
# Ici chaque pixel est moyenne sur toutes les dates/heures disponibles en 2020 (mode=:total).
matrix2 = compute_general_climatology(data_folder_ca, weight_canada, 1950, mode=:total)
# p2 = carte des moyennes temporelles par pixel en 2020.
vizumap(matrix2, weight_canada)
means_vector_calculation(matrix2, weight_canada)