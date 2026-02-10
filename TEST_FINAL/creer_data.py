import xarray as xr
import numpy as np
import pandas as pd
import os

# 8760 heures pour l'an 2000
times = pd.date_range("2000-01-01", "2000-12-31 23:00", freq="H")
lats, lons = np.linspace(41, 51, 100), np.linspace(-5, 10, 100)
temp = 15 + 10 * np.random.randn(len(times), 100, 100)

ds = xr.Dataset(
    {"t2m": (["time", "lat", "lon"], temp.astype('float32'))},
    coords={"time": times, "lat": lats, "lon": lons}
)
ds.to_netcdf("data_2000.nc")
print("Fichier data_2000.nc généré dans le dossier test.")
