import os
import time
import cdsapi
import zipfile
import xarray as xr
from pathlib import Path

# --- CONFIGURATION ---
START_YEAR = 1976
END_YEAR = 2025
OUT_DIR = Path("data/era5_land_hourly")
AREA_FR = [51.5, -5.5, 41.0, 9.8] 
GRID = [0.1, 0.1] 

# Define quarters (3 months at a time)
QUARTERS = [
    ["01", "02", "03"],
    ["04", "05", "06"],
    ["07", "08", "09"],
    ["10", "11", "12"]
]

def split_and_save(temp_file, year, months):
    """
    1. Unzips the file if needed.
    2. Opens the NetCDF.
    3. Splits it into monthly files.
    """
    work_file = temp_file # Default: assume it is already a clean NetCDF

    # --- STEP 1: Check for ZIP and Extract ---
    if zipfile.is_zipfile(temp_file):
        print(f"   📦 Detected ZIP format. Unzipping...")
        try:
            with zipfile.ZipFile(temp_file, 'r') as zip_ref:
                # We assume the zip contains one main .nc file (e.g. data.nc)
                extracted_names = zip_ref.namelist()
                if not extracted_names:
                    print("   ❌ Zip is empty.")
                    return
                
                # Extract to the same folder
                zip_ref.extractall(OUT_DIR)
                
                # The real NetCDF is the extracted file
                # Usually named 'data.nc' or similar. We find it by name.
                extracted_path = OUT_DIR / extracted_names[0]
                work_file = extracted_path
                print(f"   found extracted file: {work_file.name}")
        except Exception as e:
            print(f"   ❌ Zip extraction failed: {e}")
            return

    # --- STEP 2: Open and Split ---
    try:
        # engine='netcdf4' makes it explicit (requires pip install netCDF4)
        ds = xr.open_dataset(work_file, engine="netcdf4")
        for month in months:
            target_file = OUT_DIR / f"era5_land_fr_{year}_{month}.nc"
            
            if target_file.exists() and target_file.stat().st_size > 1000:
                 print(f"   ⏩ {month} exists. Skipping.")
                 continue

            # Select and Save
            try:
                ds_month = ds.sel(valid_time=f"{year}-{month}")
                ds_month.to_netcdf(target_file)
                print(f"   ✂️  Created: {target_file.name}")
            except KeyError:
                 print(f"   ⚠️ Month {month} not found in this chunk.")

        ds.close()

        # --- STEP 3: Cleanup ---
        # If we extracted a file (e.g. data.nc), delete it
        if work_file != temp_file:
            work_file.unlink()
        
        # Delete the original downloaded chunk (Zip or NC)
        if temp_file.exists():
            temp_file.unlink()
            
    except Exception as e:
        print(f"   ❌ Error processing NetCDF: {e}")

def main():
    c = cdsapi.Client(quiet=False)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print(f"--- 🚀 Starting Smart Batch Download (3x Faster) ---")
    
    for year in range(START_YEAR, END_YEAR + 1):
        for q_index, months in enumerate(QUARTERS):
            
            # Temporary file for the quarter
            temp_file = OUT_DIR / f"temp_quarter_{year}_Q{q_index+1}.nc"
            
            # Check if target monthly files already exist (to skip)
            if (OUT_DIR / f"era5_land_fr_{year}_{months[0]}.nc").exists():
                print(f"⏩ Skipped Quarter {q_index+1} of {year} (Files exist)")
                continue

            print(f"\n⬇️  Requesting Quarter {q_index+1} of {year} ({months})...")
            
            success = False
            while not success:
                try:
                    c.retrieve("reanalysis-era5-land", {
                        "variable": ["2m_temperature"],
                        "year": str(year),
                        "month": months, # Request 3 months at once
                        "day": [f"{d:02d}" for d in range(1, 32)],
                        "time": [f"{h:02d}:00" for h in range(24)],
                        "area": AREA_FR,
                        "grid": GRID,
                        "format": "netcdf",
                    }, str(temp_file))
                    
                    if temp_file.exists() and temp_file.stat().st_size > 1000:
                        # Success! Now split it.
                        print(f"   ✅ Downloaded Quarter. Splitting into months...")
                        split_and_save(temp_file, year, months)
                        success = True
                    else:
                        raise Exception("Downloaded file empty")
                        
                except Exception as e:
                    print(f"   ❌ Error: {e}. Retrying in 60s...")
                    time.sleep(60)

if __name__ == "__main__":
    main()