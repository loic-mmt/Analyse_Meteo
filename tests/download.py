import cdsapi
import xarray as xr
import pandas as pd
import time

# 1) Connexion à l'API CDS
c = cdsapi.Client()

# 2) Définition de la zone France (avec Corse)
#    [north, west, south, east]
area_fr = [51.5, -5.5, 41.0, 9.8]

# Paramètres de date
year = "2025"

# 3) Liste des mois à télécharger (ici: janvier et février)
months = ["01", "02"]

# 4) Liste des jours (01 à 31) – CDS gère les jours inexistants selon le mois
days = [f"{d:02d}" for d in range(1, 32)]

# 5) Liste des heures à télécharger (toutes les heures)
# -> toutes les heures
times = [f"{h:02d}:00" for h in range(24)]

# Nom du fichier NetCDF de sortie
nc_file = f"era5_t2m_{year}_{'-'.join(months)}_fr.nc"

# 4) Téléchargement ERA5 (T2m, France, plusieurs heures)
start_time = time.time()
c.retrieve(
    "reanalysis-era5-single-levels",
    {
        "product_type": "reanalysis",
        "variable": ["2m_temperature"],
        "year": year,
        "month": months,
        "day": days,
        "time": times,
        "area": area_fr,
        "grid": [0.25, 0.25],
        "format": "netcdf",
    },
    nc_file,
)
elapsed = time.time() - start_time
print(f"Téléchargement terminé en {elapsed:.1f} secondes") # 2 mois en 78 secondes

# 5) Ouverture du NetCDF et conversion en DataFrame
ds = xr.open_dataset(nc_file)


# DataFrame avec colonnes lon/lat/time/t2m
df = ds["t2m"].to_dataframe().reset_index()

# Optionnel : conversion en °C
df["t2m_c"] = df["t2m"] - 273.15


print("Terminé ! Fichier créé :")
print(" -", nc_file)


"""# 6) Export en CSV
csv_file = f"era5_t2m_{year}{month}{day}_fr_all_hours.csv"
df.to_csv(csv_file, index=False)

ds.close()

print("Terminé ! Fichiers créés :")
print(" -", nc_file)
print(" -", csv_file)
"""