import xarray as xr
import matplotlib.pyplot as plt

# 1. Chargement
ds = xr.open_dataset("data_2000.nc")

# POINT 1 : Carte fixe + Contours
plt.figure(figsize=(10, 6))
carte = ds.t2m.sel(time="2000-01-01T12:00:00", method="nearest")
carte.plot.contourf(levels=20, cmap='RdYlBu_r')
carte.plot.contour(levels=10, colors='black', alpha=0.3)
plt.title("1. Carte Température + Contours (01/01 12:00)")
plt.savefig("1_carte.png")

# POINT 2 : Moyenne annuelle
plt.figure(figsize=(10, 6))
ds.t2m.mean(dim="time").plot(cmap='viridis')
plt.title("2. Moyenne Annuelle 2000")
plt.savefig("2_moyenne.png")

# POINT 3 & 4 : Série temporelle + Poids 9x9 (Lissage)
poids_9x9 = ds.t2m.rolling(lat=9, lon=9, center=True).mean()
serie = poids_9x9.mean(dim=("lat", "lon"))
plt.figure(figsize=(12, 4))
serie.plot(color='red')
plt.title("3 & 4. Série Temporelle 2000 (Moyenne 9x9)")
plt.savefig("3_serie.png")

print("Les 3 images (1_carte.png, 2_moyenne.png, 3_serie.png) sont prêtes !")
