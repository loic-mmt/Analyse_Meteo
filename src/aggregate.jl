using NCDatasets # Pour le chargement des datasets .nc
using Glob #Pour pouvoir chercher des noms de fichiers dans un répèrtoire facillement.
using Statistics #Pour utiliser moyenne (mean) sur des matrices/dataset
using StatsPlots
using Dates
using GLM
using DataFrames
using Base.Threads

data_folder = "Analyse_Meteo/data/raw-yearly-combined/era5_fr_t2m"
data_folderRAM = "/dev/shm/era5_fr_t2m"
weight_file = "Analyse_Meteo/loic/weights_france_final.nc"
weight_bool_file = "Analyse_Meteo/src/mask_france_boolean.nc"

ds_w = NCDataset(weight_file)
ds_b = NCDataset(weight_bool_file)

weights = ds_w["final_weights"][:,:]
weights_bool = ds_b["mask"][:,:]
close(ds_w)
close(ds_b)


"""
Avec threads pour augumenter la vitesse.
"""
function visu_filtered_climatology_test(
    data_folder::String, 
    weights::Matrix{Float64}, 
    year_range; 
    selected_months::Vector{Int}=collect(1:12), 
    selected_days::Union{Integer, AbstractVector{<:Integer}, Nothing}=nothing, 
    variable_name="t2m"
)
    if selected_days isa Integer
        selected_days = [selected_days]
    end

    # Pre-transpose weights to match data layout [Lon, Lat]
    # This avoids transposing inside the loop millions of times
    weights_t = collect(weights') 

    
    println("Starting Analysis...")
    yearly_means = fill(NaN, length(year_range))
    Threads.@threads for year in year_range
        # Use a simple array for the buffer (faster than push!)
        temp_sum = 0.0
        temp_count = 0
        
        for month in selected_months
            month_str = lpad(month, 2, '0')
            files = glob("*$(year)_$(month_str)*.nc", data_folder)
            
            if isempty(files) continue end
            
            NCDataset(files[1]) do ds
                if !haskey(ds, variable_name) return end

                var = ds[variable_name]
                times = ds["valid_time"][:] 
                
                # Filter indices
                if isnothing(selected_days)
                    indices = 1:length(times)
                else
                    indices = findall(t -> day(t) in selected_days, times)
                end
                
                if isempty(indices) return end

                # Load chunk of data [Lon, Lat, Time]
                # Reading all time steps at once is usually faster for disk I/O
                data_chunk = var[:, :, indices]
                
                # --- FAST VECTORIZED CALCULATION ---
                @inbounds @simd for t in 1:size(data_chunk, 3)
                    # Get the 2D slice
                    frame = view(data_chunk, :, :, t) # 'view' avoids copying memory
                    
                    # 1. Weighted Product (Fast BLAS operation)
                    # We utilize the pre-transposed weights_t
                    w_prod = frame .* weights_t
                    
                    # 2. Filter valid data
                    # (This creates a boolean mask, but it's fast)
                    valid_mask = .!isnan.(w_prod) .& .!ismissing.(w_prod)
                    
                    if any(valid_mask)
                        # Sum only the valid parts
                        current_val = sum(w_prod[valid_mask])
                        current_weight = sum(weights_t[valid_mask])
                        
                        if current_weight > 0
                            temp_sum += (current_val / current_weight)
                            temp_count += 1
                        end
                    end
                end
            end 
        end
        
        if temp_count > 0
            mean_val = temp_sum / temp_count
            yearly_means[year-minimum(year_range)+1] = mean_val
            println("Year $year: $(round(mean_val - 273.15, digits=2))°C")
        end
    end
    # Plotting
    temps_c = yearly_means .- 273.15
    valid_years = collect(year_range)
    p = plot(valid_years, temps_c,
        title = "Filtered Climatology",
        label = "Mean Temp",
        xlabel = "Year", ylabel = "Temperature (°C)",
        lw = 2, marker = :circle, legend=:topleft
    )
    display(p)
    
    return temps_c
end
"""
Visualise le changement de température sur une période de temps donnée, au fil des années.

"""
function visu_filtered_climatology( 
    data_folder::String, 
    weights::Matrix{Float64}, 
    year_range; 
    selected_months::Vector{Int}=collect(1:12), 
    selected_days::Union{Integer, AbstractVector{<:Integer}, Nothing}=nothing, 
    variable_name="t2m"
)
    # Transformation si un seul mois 2 -> [2]
    if selected_days isa Integer
        selected_days = [selected_days]
    end
    weights = weights'
    #initilisation des mois et années.
    yearly_means = Float64[]
    valid_years = Int[]
    
    println("Starting Analysis...")
    println("Months: $selected_months")

    for year in year_range
        temp_buffer = Float64[]
        
        for month in selected_months
            month_str = lpad(month, 2, '0')
            files = glob("*$(year)_$(month_str)*.nc", data_folder)
            
            if isempty(files)
                continue
            end
            
            NCDataset(files[1]) do ds


                var = ds[variable_name]
                times = ds["valid_time"][:] 
                
                if isnothing(selected_days)
                    indices_to_keep = 1:length(times)
                else
                    indices_to_keep = findall(t -> day(t) in selected_days, times)
                end
                
                if isempty(indices_to_keep)
                    return 
                end

                data_slice = var[:, :, indices_to_keep]
                
                # On calcule la moyenne (weighted)
                for t in 1:size(data_slice)[3]
                    step_data = data_slice[:, :, t]
                    w_prod = step_data .* weights
                    
                    valid_mask = .!ismissing.(w_prod) .& .!isnan.(w_prod)
                    if any(valid_mask)
                        val = sum(w_prod[valid_mask]) / sum(weights[valid_mask])
                        push!(temp_buffer, val)
                    end
                end
            end
        end
        
        if !isempty(temp_buffer)
            push!(yearly_means, mean(temp_buffer))
            push!(valid_years, year)
            println("Year $year: $(round(mean(temp_buffer) - 273.15, digits=2))°C")
        end
    end
    extreme = extrema(valid_years)
    temps_c = yearly_means .- 273.15
    p = plot(valid_years, temps_c,
        title = "Mean temperatures in France, $extreme",
        label = "Mean Temp",
        xlabel = "Year",
        ylabel = "Temperature (°C)",
        lw = 2, marker = :circle
    )
    display(p)
    
    return temps_c
end

means = visu_filtered_climatology(data_folderRAM, weights, 1950:2025)
means = visu_filtered_climatology_test(data_folderRAM, weights, 1950:2025)

function trends_climate(means::Vector{Float64}, years_range; cutting=0)
    # Creation du dataframe pour les models et transfer en vecteur des années
    years_vec = Vector(years_range)
    df = DataFrame(Year = years_vec, Temp = means)

    # plot de base (scatter points)
    p = plot(df.Year, df.Temp,
        title = "Climate Trends Analysis",
        xlabel = "Year", ylabel = "Temperature (°C)",
        label = "Observed Mean",
        seriestype = :scatter, 
        color = :gray, alpha = 0.5,
        legend = :topleft,
        size = (800, 500)
    )

    # Calcul du model et de la prédiction
    model_global = lm(@formula(Temp ~ Year), df)
    pred_global = predict(model_global, df)

    # ajout de la tendance sur le totalité
    plot!(p, df.Year, pred_global, 
        label = "Global trend", color = :black)

    # boucle if si cutting est présent
    if cutting > minimum(years_range) && cutting < maximum(years_range)
        
        # Tri des éléments plus petits que cutting, modélisation et plot
        df1 = filter(row -> row.Year <= cutting, df)
        model1 = lm(@formula(Temp ~ Year), df1)
        pred1 = predict(model1, df1)
        plot!(p, df1.Year, pred1, label = "First trend", color = :blue)

        # Tri des éléments plus grands que cutting, modélisation et plot
        df2 = filter(row -> row.Year >= cutting, df)
        model2 = lm(@formula(Temp ~ Year), df2)
        pred2 = predict(model2, df2)
        plot!(p, df2.Year, pred2, label="Second trend", color = :red)
        
    end

    display(p)
    return p
end

trends_climate(means, 1950:2000, cutting=1980)

means = visu_filtered_climatology(data_folder, weights, 1950:2000)


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

MATRICE_MONTHS = visu_filtered_climatology_maps(data_folder, weights_bool, 1950:1990)


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

animate_climatology(MATRICE_MONTHS, 1950:1990, filename="evoltemp1950_1975.gif")


function calculate_pixel_trends(data_3d::AbstractArray{Float64, 3})
    # data_3d is (Lon, Lat, Time)
    n_lon, n_lat, n_time = size(data_3d)

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

warming = calculate_pixel_trends(MATRICE_MONTHS)

function simple_visu_trend(trend_map::AbstractArray{Float64, 2}, years::AbstractVector)

    total_change_map = trend_map .* length(years)   
    limit = maximum(abs.(filter(!isnan, total_change_map)))

    heatmap(total_change_map', 
        title = "Total Warming ($(years[1])-$(years[end]))",
        xlabel = "Longitude",
        ylabel = "Latitude",

        c = :balance,       # Blue-White-Red palette
        clims = (-limit, limit), # Forces 0 to be White
    
        aspect_ratio = :equal,
        yflip = true,       # Fix the North/South inversion
        right_margin = 5Plots.mm
    )
end

simple_visu_trend(warming, 1950:1990)

function calculate_trends_glm(data_3d::AbstractArray{Float64, 3})
    n_lon, n_lat, n_time = size(data_3d)
    
    # Initialize grids
    slope_grid = fill(NaN, n_lon, n_lat)
    p_value_grid = fill(NaN, n_lon, n_lat)
    
    time_steps = 1:n_time
    X_data = DataFrame(Time = time_steps)

    println("Calculating trends with GLM...")

    for i in 1:n_lon
        for j in 1:n_lat
            y_temps = data_3d[i, j, :]
            
            if any(isnan.(y_temps))
                continue
            end
            
            X_data.Temp = y_temps
            
            model = lm(@formula(Temp ~ Time), X_data)
            
            # Extract Results
            # coef(model)[2] is the slope (Time coefficient)
            slope_grid[i, j] = coef(model)[2]
            
            # coeftable(model).cols[4][2] is the P-value for Time
            p_value_grid[i, j] = coeftable(model).cols[4][2]
        end
    end
    
    return slope_grid, p_value_grid
end

slope, p_value = calculate_trends_glm(MATRICE_MONTHS)

function glm_visu_trend(slope_map::AbstractArray{Float64, 2}, p_map::AbstractArray{Float64, 2}, years::AbstractVector)

    # 2. Filter: Keep only significant trends (95% confidence)
    # We set non-significant pixels to NaN so they don't show up
    sig_slope_map = copy(slope_map)
    sig_slope_map[p_map .> 0.05] .= NaN

    # 3. Plot
    limit = maximum(abs.(filter(!isnan, sig_slope_map)))

    heatmap(sig_slope_map', 
        title = "Significant Warming Trends (p < 0.05)",
        c = :balance,
        clims = (-limit, limit),
        yflip = true,
        aspect_ratio = :equal
    )
end

glm_visu_trend(slope, grid, 1950:1975)


x = 1:62
y = 1:43
plot(x, y, MATRICE_MONTHS[:,:,1]',st=:surface)