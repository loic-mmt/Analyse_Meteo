using NCDatasets
using Glob
using Statistics
using Plots
using Dates

data_folder = "Analyse_Meteo/data/raw-yearly-combined/era5_fr_t2m"
weight_file = "Analyse_Meteo/loic/weights_france_final.nc"
weight_bool_file = "Analyse_Meteo/src/mask_france_boolean.nc"

ds_w = NCDataset(weight_file)
ds_b = NCDataset(weight_bool_file)
# Load weights into memory (2D array: lon x lat)
# Replace 'weights' with the actual variable name in your mask file
weights = ds_w["final_weights"][:,:]
weights_bool = ds_b["mask"][:,:]
close(ds_w)
close(ds_b)

# Calculate the sum of weights (Denominator for weighted average)
# We assume 0.0 in the weights means "outside the region"
total_weight = sum(weights)

days = Int[]
push!(days, 18)




function visu_filtered_climatology(
    data_folder::String, 
    weights::Matrix{Float64}, 
    year_range; 
    selected_months::Vector{Int}=collect(1:12), 
    selected_days::Union{Integer, AbstractVector{<:Integer}, Nothing}=nothing, 
    variable_name="t2m"
)
    # FIX 1b: Normalize inputs. If it's a single integer, turn it into a vector.
    if selected_days isa Integer
        selected_days = [selected_days]
    end

    yearly_means = Float64[]
    valid_years = Int[]
    
    println("Starting Analysis...")
    println("Months: $selected_months")
    println("Days: $(isnothing(selected_days) ? "All" : selected_days)")

    for year in year_range
        temp_buffer = Float64[]
        
        for month in selected_months
            month_str = lpad(month, 2, '0')
            files = glob("*$(year)_$(month_str)*.nc", data_folder)
            
            if isempty(files)
                continue
            end
            
            NCDataset(files[1]) do ds
                # Check if variable exists before trying to read
                if !haskey(ds, variable_name)
                    println("Warning: Var $variable_name not found in $(files[1])")
                    return # Exits the 'do' block, effectively a 'continue' for the loop
                end

                var = ds[variable_name]
                times = ds["valid_time"][:] 
                
                # --- FILTERING LOGIC ---
                if isnothing(selected_days)
                    indices_to_keep = 1:length(times)
                else
                    indices_to_keep = findall(t -> day(t) in selected_days, times)
                end
                
                # FIX 2: Use 'return' inside the 'do' block acts like 'continue' for the loop
                # If we were outside the 'do' block, we would use 'continue'.
                if isempty(indices_to_keep)
                    return 
                end

                # Extract Data
                data_slice = var[:, :, indices_to_keep]
                
                # Calculate Weighted Mean
                for t in 1:size(data_slice)[3]
                    step_data = data_slice[:, :, t]
                    w_prod = step_data .* weights' 
                    
                    valid_mask = .!ismissing.(w_prod) .& .!isnan.(w_prod)
                    if any(valid_mask)
                        val = sum(w_prod[valid_mask]) / sum(weights'[valid_mask])
                        push!(temp_buffer, val)
                    end
                end
            end # End of NCDataset do-block
        end
        
        if !isempty(temp_buffer)
            push!(yearly_means, mean(temp_buffer))
            push!(valid_years, year)
            println("Year $year: $(round(mean(temp_buffer) - 273.15, digits=2))°C")
        end
    end

    temps_c = yearly_means .- 273.15
    p = plot(valid_years, temps_c,
        title = "Filtered Climatology",
        label = "Mean Temp",
        xlabel = "Year",
        ylabel = "Temperature (°C)",
        lw = 2, marker = :circle
    )
    display(p)
    
    return temps_c
end




visu_filtered_climatology(data_folder, weights, 1950:1975, selected_months=[1,2,3])


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
                if !haskey(ds, variable_name)
                    return # Skip this file
                end

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
            clims = (min_val, max_val), # FIXED SCALING
            c = :inferno,               # Color palette
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

animate_climatology(MATRICE_MONTHS, 1950:1975, filename="evoltemp1950_1975.gif")

using Statistics

function calculate_pixel_trends(data_3d::AbstractArray{Float64, 3}, years::AbstractVector)
    # data_3d is (Lon, Lat, Time)
    n_lon, n_lat, n_time = size(data_3d)
    
    if n_time != length(years)
        error("Dimension mismatch: 3rd dim of data ($n_time) != length of years ($(length(years)))")
    end

    # Initialize the coefficient grid with NaNs
    slope_grid = fill(NaN, n_lon, n_lat)

    # Pre-calculate X variables (Time)
    # We use 1, 2, 3... instead of 1950, 1951 to keep numbers small (numerical stability)
    x = 1:n_time
    mean_x = mean(x)
    denom = sum((x .- mean_x).^2) # The bottom part of the fraction (Variance of x)

    # Iterate over every pixel
    for i in 1:n_lon
        for j in 1:n_lat
            # Extract time series for this pixel
            y = data_3d[i, j, :]
            
            # Skip if pixel is masked (contains NaNs)
            if any(isnan.(y))
                continue
            end
            
            # Calculate Slope (a)
            # a = Sum((x - mean_x) * (y - mean_y)) / Sum((x - mean_x)^2)
            mean_y = mean(y)
            numerator = sum((x .- mean_x) .* (y .- mean_y))
            
            slope = numerator / denom
            
            slope_grid[i, j] = slope
        end
    end
    
    return slope_grid
end

warming = calculate_pixel_trends(MATRICE_MONTHS, 1950:1975)
# 1. Run the calculation
# MATRICE_MONTHS is the 3D array you created in the previous step
years = 1950:1975
trend_map = calculate_pixel_trends(MATRICE_MONTHS, years)

# 2. Convert to "Total Change" (Optional but easier to read)
# Instead of "0.02°C/year", we often show "Total warming over 25 years"
total_change_map = trend_map .* length(years)

# 3. Plotting
# We calculate the max value to center the colors perfectly around 0
limit = maximum(abs.(filter(!isnan, total_change_map)))

heatmap(total_change_map', 
    title = "Total Warming ($(years[1])-$(years[end]))",
    xlabel = "Longitude",
    ylabel = "Latitude",
    
    # Color Settings
    c = :balance,       # Blue-White-Red palette
    clims = (-limit, limit), # Forces 0 to be White
    
    aspect_ratio = :equal,
    yflip = true,       # Fix the North/South inversion
    right_margin = 5Plots.mm
)