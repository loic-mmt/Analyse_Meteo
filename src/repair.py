import os
import zipfile
from pathlib import Path

# Folder containing your data
DATA_DIR = Path("data/era5_land_hourly")

def fix_files():
    files = list(DATA_DIR.glob("*.nc"))
    print(f"🔍 Found {len(files)} files. Checking for hidden ZIPs...")

    fixed_count = 0

    for file_path in files:
        # Check if it is actually a ZIP file
        if zipfile.is_zipfile(file_path):
            print(f"📦 Unzipping: {file_path.name} ...")
            
            try:
                # Open the zip
                with zipfile.ZipFile(file_path, 'r') as zip_ref:
                    # Usually contains one file like 'data.nc' or similar
                    internal_files = zip_ref.namelist()
                    if not internal_files:
                        print(f"   ⚠️ Empty zip: {file_path.name}")
                        continue
                    
                    # Extract the first file to a temporary name
                    extracted_name = internal_files[0]
                    temp_extracted_path = DATA_DIR / extracted_name
                    
                    zip_ref.extract(extracted_name, DATA_DIR)
                    
                    # Close zip so we can delete/rename
                
                # Now we have 'data.nc' (or similar) in the folder.
                # We overwrite the original 'fake .nc' (the zip) with this real .nc
                
                # 1. Delete the original zip file (which was named .nc)
                file_path.unlink()
                
                # 2. Rename the extracted file to the correct original name
                temp_extracted_path.rename(file_path)
                
                fixed_count += 1
                
            except Exception as e:
                print(f"   ❌ Error fixing {file_path.name}: {e}")
        else:
            # It's already a real NetCDF or just broken
            pass

    print(f"\n✅ Repair Complete. Fixed {fixed_count} files.")

if __name__ == "__main__":
    fix_files()




import xarray as xr




# Load the file you just unzipped manually
ds = xr.open_dataset("/mnt/data/ProjetMeteo/data/era5_land_hourly/temp_quarter_1976_Q4.nc", engine="netcdf4")


# Change these to match the file's content
year = "1976"
months = ["10", "11", "12"] 

for month in months:
    ds.sel(valid_time=f"{year}-{month}").to_netcdf(f"era5_land_fr_{year}_{month}.nc")
    print(f"Saved {month}")


