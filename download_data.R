library(reticulate)

cdsapi <- import("cdsapi")
c <- cdsapi$Client()

area_fr <- c(51.5, -5.5, 41.0, 9.8)

c$retrieve(
  "reanalysis-era5-single-levels",
  dict(
    product_type = "reanalysis",
    variable     = list("2m_temperature"),
    year         = "2025",
    month        = "01",
    day          = "01",
    time         = "00:00",
    area         = list(area_fr[1], area_fr[2], area_fr[3], area_fr[4]),
    grid         = list(0.25, 0.25),
    format       = "netcdf"
  ),
  "test_era5_t2m_20250101_0000.nc"
)

library(terra)
library(sf)
library(geodata)

nc_file <- "test_era5_t2m_20250101_0000.nc"

# Charge le raster ERA5 (lon/lat, une seule heure ici)
r <- rast(nc_file)

# Récupère les régions françaises (niveau 1) depuis GADM
fr_admin1 <- geodata::gadm("FRA", level = 1, path = "data/gadm")
fr_admin1_sf <- st_as_sf(fr_admin1)

# On enlève la Corse
fr_mainland_sf <- subset(fr_admin1_sf, !NAME_1 %in% c("Corse"))
fr_mainland <- vect(fr_mainland_sf)  # SpatVector pour terra

# Découpe au contour de la France métropolitaine SANS Corse
r_crop <- crop(r, fr_mainland)
r_mask <- mask(r_crop, fr_mainland)

# Sauvegarde éventuelle en NetCDF propre (optionnel)
writeRaster(r_mask, "test_era5_t2m_20250101_0000_france_sans_corse.nc",
            overwrite = TRUE)

df <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)

# Renomme proprement (le nom exact de la colonne dépendra du fichier)
colnames(df) <- c("lon", "lat", "t2m")

write.csv(df,
          "test_era5_t2m_20250101_0000_france_sans_corse.csv",
          row.names = FALSE)