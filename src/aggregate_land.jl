using NCDatasets # Pour le chargement des datasets .nc
using Glob #Pour pouvoir chercher des noms de fichiers dans un répèrtoire facillement.
using Statistics #Pour utiliser moyenne (mean) sur des matrices/dataset
using StatsPlots
using Dates
using GLM
using DataFrames
using Base.Threads
using RollingFunctions

data_folder = "/mnt/data/ProjetMeteo/data/era5_land_hourly"
data_folderRAM = "/dev/shm/era5_fr_t2m"
#weight_file = "Analyse_Meteo/loic/weights_france_final.nc"
#weight_bool_file = "Analyse_Meteo/src/mask_france_boolean.nc"
weight_file ="/mnt/data/ProjetMeteo/Analyse_Meteo/loic/weights_france_final.nc"
weight_bool_file = "/mnt/data/ProjetMeteo/Analyse_Meteo/src/mask_france_boolean_land.nc"
ds_w = NCDataset(weight_file)
ds_b = NCDataset(weight_bool_file)

weights = ds_w["final_weights"][:,:]
weights_bool = ds_b["mask"][:,:]
close(ds_w)
close(ds_b)

function visu_filtered_climatology_maps(
    data_folder::String, 
    weights::Matrix{Float64}, 
    year_range; 
    selected_months::Vector{Int}=collect(1:12), 
    selected_days::Union{Integer, AbstractVector{<:Integer}, Nothing}=nothing, 
    variable_name="t2m"
)
    # 1. Normalize inputs
    if selected_days isa Integer
        selected_days = [selected_days]
    end

    visual_mask = fill(NaN, size(weights)) 
    visual_mask[weights .> 0.9] .= 1.0

    # We store the 2D maps for each year in a list first
    yearly_maps_list = Matrix{Float64}[]
    valid_years = Int[]
    
    println("Starting Map Analysis...")

    for year in year_range
        # Temporary grids for the current year
        year_sum_grid = nothing
        year_count_grid = nothing
        
        for month in selected_months
            month_str = lpad(month, 2, '0')
            files = glob("*$(year)_$(month_str)*.nc", data_folder)
            
            if isempty(files)
                continue
            end
            
            NCDataset(files[1]) do ds

                var = ds[variable_name]
                times = ds["valid_time"][:] 
                
                # Filter days
                if isnothing(selected_days)
                    indices_to_keep = 1:length(times)
                else
                    indices_to_keep = findall(t -> day(t) in selected_days, times)
                end
                
                if isempty(indices_to_keep)
                    return # Skip this file
                end

                # Extract Data Slice [Lon, Lat, Time]
                data_slice = var[:, :, indices_to_keep]
                
                # Initialize grids on the first valid file of the year
                if isnothing(year_sum_grid)
                    dims = size(data_slice)
                    year_sum_grid = zeros(Float64, dims[1], dims[2])
                    year_count_grid = zeros(Int, dims[1], dims[2])
                end

                # Accumulate
                for t in 1:size(data_slice)[3]
                    frame = data_slice[:, :, t]
                    
                    # Only sum valid numbers (skip NaN/Missing)
                    valid_mask = .!ismissing.(frame) .& .!isnan.(frame)
                    
                    if any(valid_mask)
                        year_sum_grid[valid_mask] .+= frame[valid_mask]
                        year_count_grid[valid_mask] .+= 1
                    end
                end
            end 
        end
        
        # End of Year: Compute Mean
        if !isnothing(year_sum_grid) && any(year_count_grid .> 0)
            mean_grid = fill(NaN, size(year_sum_grid))
            valid_pixels = year_count_grid .> 0
            
            # Sum / Count
            mean_grid[valid_pixels] .= year_sum_grid[valid_pixels] ./ year_count_grid[valid_pixels]
            mean_grid = mean_grid .* visual_mask'
            # Convert to Celsius
            mean_grid .-= 273.15
            
            push!(yearly_maps_list, mean_grid)
            push!(valid_years, year)
            println("Year $year processed.")
        end
    end

    # 2. Convert Vector of Matrices -> 3D Array [Lon, Lat, Year]
    if isempty(yearly_maps_list)
        println("No data found!")
        return Array{Float64}(undef, 0, 0, 0)
    end

    # 'cat' stacks the matrices along dimension 3
    final_3d_matrix = cat(yearly_maps_list..., dims=3)
    
    println("Output dimensions: $(size(final_3d_matrix)) (Lon x Lat x Years)")
    
    return final_3d_matrix
end

MATRICE_MONTHS = visu_filtered_climatology_maps(data_folder, weights_bool, 1950:1975)

function animate_climatology(data_3d::AbstractArray{Float64, 3}, valid_years::AbstractVector; filename="temperature_evolution.gif")
    
    println("Generating animation...")

    # 1. Determine fixed color limits for the whole period
    # We ignore NaNs so they don't break the min/max calculation
    valid_data = filter(!isnan, data_3d)
    if isempty(valid_data)
        println("Error: Data contains only NaNs.")
        return
    end
    min_val, max_val = minimum(valid_data), maximum(valid_data)
    
    # 2. Create the Animation object
    anim = @animate for i in 1:length(valid_years)
        year = valid_years[i]
        
        # Extract the 2D map for this year
        # Transpose (') is usually needed because Julia arrays are Col-Major
        # but heatmap expects [x, y]. 
        current_map = data_3d[:, :, i]'
        
        heatmap(current_map,
            title = "Mean Temperature: $year",
            clims = (min_val, max_val),
            c = :thermal,   # Color palette
            xlabel = "Longitude",
            ylabel = "Latitude",
            aspect_ratio = :equal,
            right_margin = 5Plots.mm,
            yflip = true    # Give space for the colorbar
        )
    end

    # 3. Save the GIF
    # fps = frames per second. 
    gif(anim, filename, fps = 5) 
    println("Saved animation to $filename")
end

animate_climatology(MATRICE_MONTHS, 1950:1975, filename="evoltemp1950_2025thermal.gif")


function visu_filtered_climatology_maps(
    data_folder::String, 
    weights::Matrix{Float64}, 
    year_range; 
    selected_months::Vector{Int}=collect(1:12), 
    selected_days::Union{Integer, AbstractVector{<:Integer}, Nothing}=nothing, 
    variable_name="t2m",
    export_path::Union{String, Nothing}=nothing
)
    # 1. Normalize inputs
    if selected_days isa Integer
        selected_days = [selected_days]
    end

    visual_mask = fill(NaN, size(weights)) 
    visual_mask[weights .> 0.9] .= 1.0

    yearly_maps_list = Matrix{Float64}[]
    valid_years = Int[]
    
    # We need to capture coordinates for the export
    lons_vector = nothing
    lats_vector = nothing
    
    println("Starting Map Analysis...")

    for year in year_range
        year_sum_grid = nothing
        year_count_grid = nothing
        
        for month in selected_months
            month_str = lpad(month, 2, '0')
            files = glob("*$(year)_$(month_str)*.nc", data_folder)
            
            if isempty(files) continue end
            
            NCDataset(files[1]) do ds
                # --- CAPTURE COORDINATES (ONCE) ---
                if isnothing(lons_vector)
                    lons_vector = ds["longitude"][:]
                    lats_vector = ds["latitude"][:]
                end

                var = ds[variable_name]
                times = ds["valid_time"][:] 
                
                # Filter days
                if isnothing(selected_days)
                    indices_to_keep = 1:length(times)
                else
                    indices_to_keep = findall(t -> day(t) in selected_days, times)
                end
                
                if isempty(indices_to_keep) return end

                data_slice = var[:, :, indices_to_keep]
                
                if isnothing(year_sum_grid)
                    dims = size(data_slice)
                    year_sum_grid = zeros(Float64, dims[1], dims[2])
                    year_count_grid = zeros(Int, dims[1], dims[2])
                end

                # Accumulate
                for t in 1:size(data_slice)[3]
                    frame = data_slice[:, :, t]
                    valid_mask = .!ismissing.(frame) .& .!isnan.(frame)
                    
                    if any(valid_mask)
                        year_sum_grid[valid_mask] .+= frame[valid_mask]
                        year_count_grid[valid_mask] .+= 1
                    end
                end
            end 
        end
        
        # End of Year Processing
        if !isnothing(year_sum_grid) && any(year_count_grid .> 0)
            mean_grid = fill(NaN, size(year_sum_grid))
            valid_pixels = year_count_grid .> 0
            
            mean_grid[valid_pixels] .= year_sum_grid[valid_pixels] ./ year_count_grid[valid_pixels]
            mean_grid = mean_grid .* visual_mask' # Apply mask
            mean_grid .-= 273.15 # Kelvin -> Celsius
            
            push!(yearly_maps_list, mean_grid)
            push!(valid_years, year)
            println("Year $year processed.")
        end
    end

    if isempty(yearly_maps_list)
        println("No data found!")
        return Array{Float64}(undef, 0, 0, 0)
    end

    # Stack to 3D Matrix [Lon, Lat, Time]
    final_3d_matrix = cat(yearly_maps_list..., dims=3)
    println("Final dimensions: $(size(final_3d_matrix))")

    # --- 2. EXPORT LOGIC (NetCDF) ---
    if !isnothing(export_path)
        println("Exporting to $export_path ...")
        
        NCDataset(export_path, "c") do ds_out
            # Define Dimensions
            defDim(ds_out, "longitude", length(lons_vector))
            defDim(ds_out, "latitude", length(lats_vector))
            defDim(ds_out, "year", length(valid_years))
            
            # Define Variables (Coordinates)
            v_lon = defVar(ds_out, "longitude", Float64, ("longitude",))
            v_lon[:] = lons_vector
            
            v_lat = defVar(ds_out, "latitude", Float64, ("latitude",))
            v_lat[:] = lats_vector
            
            v_year = defVar(ds_out, "year", Int32, ("year",))
            v_year[:] = valid_years
            
            # --- THE FIX IS HERE ---
            # We define _FillValue and units INSIDE defVar, before writing any data
            v_temp = defVar(ds_out, "temperature", Float64, ("longitude", "latitude", "year"), attrib = [
                "_FillValue" => NaN,
                "units" => "Celsius",
                "long_name" => "Mean Temperature"
            ])
            
            # Write data AFTER defining attributes
            v_temp[:,:,:] = final_3d_matrix
            
            # Global attributes can be set anytime
            ds_out.attrib["title"] = "Processed Climatology"
        end
        println("Export successful")
    end
    
    return final_3d_matrix
end


matrix = visu_filtered_climatology_maps(
    data_folder, 
    weights_bool, 
    1950:2025,
    export_path="output_climate_land.nc"
)