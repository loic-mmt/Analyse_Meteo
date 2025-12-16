import os
import time
from pathlib import Path

import cdsapi
import xarray as xr


# Chunked by month to avoid huge single requests/files


# 1) Connexion à l'API CDS
c = cdsapi.Client()

# 2) Définition de la zone France (avec Corse) [north, west, south, east]

AREA_FR = [51.5, -5.5, 41.0, 9.8]
GRID = [0.25, 0.25]

# 3) Période à télécharger (1950-2024 inclus)
START_YEAR = 2013
END_YEAR = 2025
YEARS = list(range(START_YEAR, END_YEAR + 1))
MONTHS = list(range(1, 13))
DAYS = [f"{d:02d}" for d in range(1, 32)]  # CDS gère les jours inexistants selon le mois
TIMES = [f"{h:02d}:00" for h in range(24)]

# 4) Dossier de sortie
OUT_DIR = Path("era5_fr_t2m")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def retrieve_one_month(year: int, month: int, out_path: Path) -> None:
    """Download one month of ERA5 t2m for France to NetCDF."""
    request = {
        "product_type": "reanalysis",
        "variable": ["2m_temperature"],
        "year": f"{year:04d}",
        "month": f"{month:02d}",
        "day": DAYS,
        "time": TIMES,
        "area": AREA_FR,
        "grid": GRID,
        "format": "netcdf",
    }

    start = time.time()
    c.retrieve(
        "reanalysis-era5-single-levels",
        request,
        str(out_path),
    )
    elapsed = time.time() - start
    print(f"OK {year}-{month:02d} -> {out_path.name} ({elapsed:.1f}s)")


def main() -> None:
    total_start = time.time()
    downloaded = 0
    skipped = 0

    for year in YEARS:
        for month in MONTHS:
            out_file = OUT_DIR / f"era5_t2m_fr_{year:04d}_{month:02d}.nc"

            # Skip if already present
            if out_file.exists() and out_file.stat().st_size > 0:
                skipped += 1
                continue

            # Retry a couple of times (network/CDS hiccups happen)
            for attempt in range(1, 4):
                try:
                    retrieve_one_month(year, month, out_file)
                    downloaded += 1
                    break
                except Exception as e:
                    if attempt == 3:
                        raise
                    print(f"WARN {year}-{month:02d} attempt {attempt} failed: {e}")
                    time.sleep(5)

    total_elapsed = time.time() - total_start
    print("\nTerminé !")
    print(f"Dossier : {OUT_DIR.resolve()}")
    print(f"Téléchargés : {downloaded}")
    print(f"Déjà présents : {skipped}")
    print(f"Temps total : {total_elapsed/60:.1f} min")

    # IMPORTANT:
    # Ne pas faire ds = xr.open_dataset(...) sur un fichier unique énorme.
    # Pour analyser ensuite, utiliser:
    #   ds = xr.open_mfdataset(str(OUT_DIR / 'era5_t2m_fr_*.nc'), combine='by_coords')


if __name__ == "__main__":
    main()