import cdsapi
import xarray as xr
import pandas as pd

# 1) Connexion à l'API CDS
c = cdsapi.Client()

# 2) Définition de la zone France (avec Corse)
#    [north, west, south, east]
area_fr = [51.5, -5.5, 41.0, 9.8]

# Paramètres de date
year = "2025"
month = "01"
day = "01"

# 3) Liste des heures à télécharger
# -> toutes les heures
times = [f"{h:02d}:00" for h in range(24)]
# times = ["00:00", "01:00", "02:00", "03:00"]

# Nom du fichier NetCDF de sortie
nc_file = f"era5_t2m_{year}{month}{day}_fr.nc"

# 4) Téléchargement ERA5 (T2m, France, plusieurs heures)
c.retrieve(
    "reanalysis-era5-single-levels",
    {
        "product_type": "reanalysis",
        "variable": ["2m_temperature"],
        "year": year,
        "month": month,
        "day": day,
        "time": times,
        "area": area_fr,
        "grid": [0.25, 0.25],
        "format": "netcdf",
    },
    nc_file,
)

# 5) Ouverture du NetCDF et conversion en DataFrame
ds = xr.open_dataset(nc_file)


# DataFrame avec colonnes lon/lat/time/t2m
df = ds["t2m"].to_dataframe().reset_index()

# Optionnel : conversion en °C
df["t2m_c"] = df["t2m"] - 273.15

# 6) Export en CSV
csv_file = f"era5_t2m_{year}{month}{day}_fr_all_hours.csv"
df.to_csv(csv_file, index=False)

ds.close()

print("Terminé ! Fichiers créés :")
print(" -", nc_file)
print(" -", csv_file)