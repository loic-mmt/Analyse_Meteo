using JLD2
using Glob
using NCDatasets
using Dates
using Plots
using Base.Threads
using Statistics


function convert_to_fast_format(data_folder, output_file)

    files = glob("*.nc", data_folder)

    # 1. Capture Coordinates from the first file
    first_ds = NCDataset(files[1])
    lats = first_ds["latitude"][:]
    lons = first_ds["longitude"][:]
    close(first_ds)
    
    # 2. Sort files chronologically (to ensure correct order)
    # Assumes filename format like "era5_1950_01.nc"
    sort!(files, by = f -> begin
        m = match(r"_(\d{4})_(\d{2})", basename(f))
        parse(Int, m[1]) * 100 + parse(Int, m[2]) 
    end)
    
    # 3. Prepare Containers
    # We store a VECTOR of matrices (one chunk per original file)
    all_temps = Vector{Array{Float32, 3}}() 
    all_time_vectors = Vector{Vector{Any}}() # Stores the exact time axis for each chunk
    all_years = Int[]
    all_months = Int[]
    
    println("Processing $(length(files)) files...")
    
    for (i, f) in enumerate(files)
        NCDataset(f, "r") do ds
            # Extract Year/Month from filename for quick indexing later
            m = match(r"_(\d{4})_(\d{2})", basename(f))
            year, month = parse(Int, m[1]), parse(Int, m[2])
            
            # Check for variable (usually "t2m" or "t2m_0001")
            # You can adapt "t2m" if your variable has a different name
            var_name = haskey(ds, "t2m") ? "t2m" : first(keys(ds)) 
            
            if haskey(ds, var_name)
                # Read Data (Load into RAM)
                var = ds[var_name][:,:,:]
                
                # Read Time Axis (Critical for day filtering)
                times = ds["valid_time"][:] 

                # Store converted to Float32 (saves 50% RAM/Disk compared to Float64)
                # TIP: Use Float16.(var) here if you want to reduce the file size to ~3.5GB
                push!(all_temps, Float16.(var))
                
                push!(all_time_vectors, times)
                push!(all_years, year)
                push!(all_months, month)
            end
        end
        if i % 10 == 0 print(".") end
    end
    
    
    # Save everything using LZ4 compression (Fast & Efficient)
    jldsave(output_file; 
        lats, lons, 
        all_temps, 
        all_time_vectors, 
        all_years, 
        all_months, 
        compress=true)
    
end



# Inside compute_general_climatology...

if endswith(data_folder, ".jld2")
    println("⚡ Detected JLD2 Fast Format.")
    (sums, counts, lons, lats) = accumulate_data_jld2(
        data_folder, year_range, mode, selected_months, selected_days, variable_name
    )
else
    # Fallback to NetCDF folder logic
    (sums, counts, lons, lats) = accumulate_data(
        data_folder, year_range, mode, selected_months, selected_days, variable_name
    )
end

function accumulate_data_jld2(jld2_path::String, year_range, mode, months, days, var_name)
    sums_dict = Dict{Any, Matrix{Float64}}()
    counts_dict = Dict{Any, Matrix{Int}}()
    lons, lats = nothing, nothing

    println("Loading JLD2 file structure...")
    
    jldopen(jld2_path, "r") do file
        # 1. Load Metadata (Lightweight, instant)
        lons = file["lons"]
        lats = file["lats"]
        all_years = file["all_years"]
        all_months = file["all_months"]
        
        n_chunks = length(all_years)
        println("File contains $n_chunks chunks (months). Processing...")

        # 2. Iterate through chunks
        for i in 1:n_chunks
            current_year = all_years[i]
            current_month = all_months[i]

            # OPTIMIZATION: Skip chunks that are not in the requested range
            # This makes the function extremely fast for subsetting
            if !(current_year in year_range) || !(current_month in months)
                continue
            end

            # 3. Load ONLY this chunk into RAM
            # file["key"][i] reads only the i-th element of the vector from disk
            chunk_data = file["all_temps"][i]       
            chunk_times = file["all_time_vectors"][i] 

            # 4. Filter Days (Using the time vector we saved)
            if isnothing(days)
                indices = 1:length(chunk_times)
            else
                indices = findall(t -> day(t) in days, chunk_times)
            end

            # 5. Process the valid slices
            if !isempty(indices)
                data_slice = chunk_data[:, :, indices]
                
                for t in 1:length(indices)
                    time_val = chunk_times[indices[t]]
                    
                    # Determine Aggregation Key (Monthly, Yearly, or Total)
                    key = get_binning_key(mode, current_year, current_month, time_val)

                    # Init Dict if new key
                    if !haskey(sums_dict, key)
                        dims = size(data_slice)[1:2]
                        sums_dict[key] = zeros(Float64, dims)
                        counts_dict[key] = zeros(Int, dims)
                    end

                    # Accumulate
                    frame = data_slice[:, :, t]
                    
                    # Check validity (NaN check covers both missing and NaN for Float32)
                    valid = .!isnan.(frame) 
                    
                    if any(valid)
                        sums_dict[key][valid] .+= frame[valid]
                        counts_dict[key][valid] .+= 1
                    end
                end
            end
            # Optional: nice progress bar logic could go here
        end
    end
    
    println("JLD2 processing complete.")
    return sums_dict, counts_dict, lons, lats
end