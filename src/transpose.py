import xarray as xr
import glob
import os

def batch_fix_dimensions(folder_path):
    # 1. Find all .nc files in the folder
    # construct the pattern (e.g., "data/*.nc")
    search_pattern = os.path.join(folder_path, "*.nc")
    files = glob.glob(search_pattern)
    
    print(f"📂 Found {len(files)} files in '{folder_path}'")

    for file_path in files:
        try:
            print(f"   🔄 Processing: {os.path.basename(file_path)}...", end="")
            
            # 2. Open the file
            with xr.open_dataset(file_path) as ds:
                
                # 3. Transpose (Reorder Dimensions)
                # Force order: Time -> Latitude -> Longitude
                # (Missing dimensions are ignored, existing ones are reordered)
                ds_fixed = ds.transpose("time", "latitude", "longitude", ...)
                
                # 4. Save to a Temporary File
                # We CANNOT overwrite the file while it is open, so we save to a temp name.
                temp_path = file_path.replace(".nc", "_temp.nc")
                ds_fixed.to_netcdf(temp_path)
            
            # 5. Overwrite the Original
            # Now that 'ds' is closed (thanks to 'with'), we can replace the file.
            os.replace(temp_path, file_path)
            
            print(" ✅ Done.")

        except Exception as e:
            print(f" ❌ Error: {e}")
            # Clean up temp file if it was created but failed later
            if 'temp_path' in locals() and os.path.exists(temp_path):
                os.remove(temp_path)

# --- CONFIGURATION ---
# Change this to your specific folder path
TARGET_FOLDER = "data/masks" 

if __name__ == "__main__":
    batch_fix_dimensions(TARGET_FOLDER)