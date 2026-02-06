import xarray as xr
import glob
import os
import sys

def concatenate_to_global_file(source_folder, output_file, year_range=None):
    """
    Concatenates monthly NetCDF files into a single global file using Xarray,
    optimized for fast random access reading in Julia/Python.
    """
    print(f"🚀 Searching for files in: {source_folder}")
    
    # 1. Get List of Files
    files = []
    if year_range:
        for year in year_range:
            pattern = os.path.join(source_folder, f"*{year}_*.nc")
            found = glob.glob(pattern)
            files.extend(found)
    else:
        files = glob.glob(os.path.join(source_folder, "*.nc"))
    
    files.sort()
    
    if not files:
        print("❌ Error: No files found! Check your path or year range.")
        return

    print(f"📦 Found {len(files)} files to merge.")
    print("⏳ Opening datasets...")

    # 2. Open Files with "Auto" Chunking
    # We do NOT use -1 here. We let Dask manage the READ memory automatically.
    try:
        ds = xr.open_mfdataset(
            files, 
            chunks={"time": "auto"},  # Let Dask decide the safe READ size
            combine="by_coords",
            parallel=True,
            engine="h5netcdf"
        )
    except Exception as e:
        print(f"❌ Error opening files: {e}")
        return

    print(f"📊 Dataset Dimensions: {ds.dims}")
    print(f"⚙️  Configuring optimized chunking for Disk...")

    # 3. Define OUTPUT Encoding (The Speed Secret)
    # We force the file to be saved in small "bricks" (Chunks) on the hard drive.
    # This allows Julia to pull out just 1 brick (e.g., 1 week) without reading the whole wall.
    
    # Calculate spatial size
    n_lat = ds.dims['latitude']
    n_lon = ds.dims['longitude']
    
    # We set time chunk to ~100 steps. 
    # For hourly data, this is ~4 days per block. Very fast for random access.
    # For monthly data, this is ~8 years per block.
    time_chunk_size = 100 

    encoding = {}
    for var in ds.data_vars:
        encoding[var] = {
            "zlib": True,        # Enable Compression
            "complevel": 1,      # Level 1 is FAST (Level 9 is slow)
            "shuffle": True,     # Improves compression speed
            # (Time, Lat, Lon) -> Writes full maps, but cuts time into slices
            "chunksizes": (time_chunk_size, n_lat, n_lon) 
        }

    print(f"💾 Saving to {output_file} (Optimized)...")

    # 4. Write to NetCDF
    try:
        ds.to_netcdf(
            output_file, 
            compute=True,
            format="NETCDF4",
            engine="h5netcdf",
            encoding=encoding  # <--- Apply the optimization
        )
        print(f"✅ Success! File saved: {output_file}")
    
    except MemoryError:
        print("❌ RAM Error: The system ran out of memory.")
    except Exception as e:
        print(f"❌ Save Failed: {e}")

# --- CONFIGURATION ---
if __name__ == "__main__":
    # ⚠️ UPDATE YOUR PATHS HERE
    # Use the path where your monthly files are located
    SOURCE_PATH = "/mnt/data/ProjetMeteo/Analyse_Meteo/data/raw_monthly_combined/basic" 
    OUTPUT_PATH = "era5_global_1950_2025_basic.nc"

    # Run the function
    concatenate_to_global_file(SOURCE_PATH, OUTPUT_PATH, range(1950, 2026))