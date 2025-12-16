using JLD2
using Glob
using NCDatasets
using Dates
using Plots
using Base.Threads
using Statistics

function convert_to_fast_format(data_folder, output_file)
    println("--- Starting Conversion to JLD2 ---")
    
    files = glob("*.nc", data_folder)
    
    # We will store everything in big arrays
    # Assuming all files have the same Lat/Lon dimensions
    first_ds = NCDataset(files[1])
    lats = first_ds["latitude"][:]
    lons = first_ds["longitude"][:]
    close(first_ds)
    
    n_lat, n_lon = length(lats), length(lons)
    
    # Sort files by date to ensure time order
    # Extract date from filename, e.g., "era5_1950_01.nc"
    sort!(files, by = f -> begin
        m = match(r"_(\d{4})_(\d{2})", basename(f))
        parse(Int, m[1]) * 100 + parse(Int, m[2]) # e.g. 195001
    end)
    
    # Prepare data containers
    # We use Float32 (4 bytes) instead of Float64 (8 bytes) to save 50% disk/RAM
    all_temps = Vector{Array{Float32, 3}}() 
    all_years = Int[]
    all_months = Int[]
    
    println("Processing $(length(files)) files...")
    
    for (i, f) in enumerate(files)
        NCDataset(f, "r") do ds
            # Extract Year/Month from filename or metadata
            m = match(r"_(\d{4})_(\d{2})", basename(f))
            year, month = parse(Int, m[1]), parse(Int, m[2])
            
            # Read Variable (Load into RAM immediately)
            # Adjust "t2m" if your variable name is different
            if haskey(ds, "t2m")
                var = ds["t2m"][:,:,:] 
                
                # Convert to Float32 immediately to reduce size
                push!(all_temps, Float32.(var))
                push!(all_years, year)
                push!(all_months, month)
            end
        end
        if i % 10 == 0 print(".") end
    end
    
    println("\nSaving to $output_file ... (This may take a minute)")
    
    # Save to JLD2
    # compress=true uses LZ4, which is insanely fast to decompress compared to NetCDF
    jldsave(output_file; lats, lons, all_temps, all_years, all_months, compress=true)
    
    println("\n✅ Done! File saved: $output_file")
end

# Run it
convert_to_fast_format("Analyse_Meteo/data/raw-yearly-combined/era5_fr_t2m", "/home/arthur/Desktop/DATAJLD2.jld2")


weight_file = "Analyse_Meteo/loic/weights_france_final.nc"
ds_w = NCDataset(weight_file)
weights = ds_w["final_weights"][:,:]
close(ds_w)

function visu_jld2(jld_file, weights::Matrix{Float64}, year_range)
    println("Loading JLD2 Data...")
    
    # 1. Load Everything (Super Fast)
    # The data is essentially "dumped" from disk to RAM directly
    data = load(jld_file)
    
    all_temps = data["all_temps"]   # Vector of 3D Arrays
    all_years = data["all_years"]
    
    # Transpose weights once
    weights_t = collect(weights')
    
    # 2. Filter Years
    years_vec = collect(year_range)
    n_years = length(years_vec)
    yearly_means = fill(NaN, n_years)
    
    println("Computing on $(Threads.nthreads()) threads...")
    
    # 3. Parallel Compute
    Threads.@threads for i in 1:n_years
        target_year = years_vec[i]
        
        # Find all indices in the data that match this year
        # (This is fast because all_years is just a list of integers)
        indices = findall(y -> y == target_year, all_years)
        
        if isempty(indices) continue end
        
        temp_sum = 0.0
        temp_count = 0
        
        # Iterate over the months (chunks) for this year
        for idx in indices
            chunk = all_temps[idx] # Get the 3D array [Lon, Lat, Time]
            
            # Inner Loop (Vectorized)
            for t in 1:size(chunk, 3)
                frame = view(chunk, :, :, t)
                
                w_prod = frame .* weights_t
                # Float32 can have NaNs too
                valid_mask = .!isnan.(w_prod)
                
                if any(valid_mask)
                    # Compute in Float64 for precision
                    val = sum(Float64, w_prod[valid_mask])
                    w = sum(Float64, weights_t[valid_mask])
                    
                    if w > 0
                        temp_sum += (val / w)
                        temp_count += 1
                    end
                end
            end
        end
        
        if temp_count > 0
            yearly_means[i] = temp_sum / temp_count
        end
    end
    
    # 4. Plot
    mask = .!isnan.(yearly_means)
    final_years = years_vec[mask]
    final_temps = yearly_means[mask] .- 273.15
    
    p = plot(final_years, final_temps, marker=:circle, label="Temp JLD2")
    display(p)
    return(final_temps)
end

visu_jld2("/home/arthur/Desktop/DATAJLD2.jld2", weights, 1950:2025)

using JLD2
using Plots
using Base.Threads
using Statistics

function visu_fastest(jld_file, weights::Matrix{Float64}, year_range)
    # 1. Load Data (Instant)
    println("Loading JLD2...")
    data = load(jld_file)
    all_temps = data["all_temps"]
    all_years = data["all_years"]
    
    # Flatten weights to 1D array for faster iteration
    # (Computer memory is linear, so 1D is always faster than 2D)
    weights_flat = vec(collect(weights')) 
    
    years_vec = collect(year_range)
    n_years = length(years_vec)
    yearly_means = fill(NaN, n_years)

    println("Computing on $(Threads.nthreads()) threads (Zero-Alloc Mode)...")

    # 2. Parallel Loop
    Threads.@threads for i in 1:n_years
        target_year = years_vec[i]
        
        # Fast lookup using searchsorted (if sorted) or just findall
        indices = findall(y -> y == target_year, all_years)
        if isempty(indices) continue end
        
        temp_sum = 0.0
        temp_count = 0
        
        for idx in indices
            chunk = all_temps[idx] # [Lon, Lat, Time]
            
            # We iterate strictly over the Time dimension
            n_times = size(chunk, 3)
            
            for t in 1:n_times
                # Access the 2D frame as a linear 1D slice (Zero cost)
                # We interpret the memory as a simple list of numbers
                frame_start = (t-1) * length(weights_flat) + 1
                frame_end = t * length(weights_flat)
                
                # Reshape/View is not even needed if we use linear indexing on the parent array
                # But 'view' is safer and usually optimized away
                frame_flat = view(chunk, :, :, t)
                
                # --- THE KERNEL (Optimized) ---
                numerator = 0.0
                denominator = 0.0
                
                # @inbounds: Don't check if index exists (Dangerous but fast)
                # @simd: Calculate 8 numbers at once
                @inbounds @simd for k in eachindex(weights_flat)
                    val = frame_flat[k]
                    w = weights_flat[k]
                    
                    # Check NaN on the value (assuming weights are always valid)
                    if !isnan(val)
                        numerator += val * w
                        denominator += w
                    end
                end
                
                if denominator > 0
                    temp_sum += (numerator / denominator)
                    temp_count += 1
                end
            end
        end
        
        if temp_count > 0
            yearly_means[i] = temp_sum / temp_count
        end
    end
    
    # 3. Plot
    mask = .!isnan.(yearly_means)
    final_years = years_vec[mask]
    final_temps = yearly_means[mask] .- 273.15
    
    p = plot(final_years, final_temps, 
        marker=:circle, 
        label="Temp (Optimized)", 
        title="Climatology (Zero Alloc)",
        lw=2
    )
    display(p)
end