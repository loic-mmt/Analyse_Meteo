using NCDatasets
using Glob
using Statistics
using Plots
using Dates

data_folder = "Analyse_Meteo/data/raw-yearly-combined/era5_fr_t2m"
weight_file = "Analyse_Meteo/loic/weights_france_final.nc"

ds_w = NCDataset(weight_file)
# Load weights into memory (2D array: lon x lat)
# Replace 'weights' with the actual variable name in your mask file
weights = ds_w["final_weights"][:,:] 
close(ds_w)

# Calculate the sum of weights (Denominator for weighted average)
# We assume 0.0 in the weights means "outside the region"
total_weight = sum(weights)




function visu_filtered_climatology(
    data_folder::String, 
    weights::Matrix{Float64}, 
    year_range; 
    selected_months::Vector{Int}=collect(1:12), # Default: All months
    selected_days::Union{AbstractVector{Int}, Nothing}=nothing, # Default: All days
    variable_name="t2m"
)
    # 1. Prepare Storage
    # We store one value per year (the average of the selected period for that year)
    yearly_means = Float64[]
    valid_years = Int[]
    

    println("Starting Analysis...")
    println("Months: $selected_months")
    println("Days: $(isnothing(selected_days) ? "All" : selected_days)")

    for year in year_range
        temp_buffer = Float64[]
        
        # Loop only over the months you asked for
        for month in selected_months
            month_str = lpad(month, 2, '0')
            files = glob("*$(year)_$(month_str)*.nc", data_folder)
            
            if isempty(files)
                continue
            end
            
            # Open file
            NCDataset(files[1]) do ds
                var = ds[variable_name]
                times = ds["valid_time"][:] # Load time vector
                
                # --- FILTERING LOGIC ---
                # 1. Identify indices to keep
                if isnothing(selected_days)
                    # Keep all time steps in this file
                    indices_to_keep = 1:length(times)
                else
                    # Find indices where the day matches your request
                    # e.g., only days 1, 2, ..., 10
                    indices_to_keep = findall(t -> day(t) in selected_days, times)
                end
                
                # If no days match (e.g., looking for day 31 in Feb), skip
                if isempty(indices_to_keep)
                    return
                end

                # 2. Extract Data (Only load what we need)
                # Syntax: var[:, :, indices]
                data_slice = var[:, :, indices_to_keep]
                
                # 3. Calculate Weighted Mean for this slice
                # We iterate over the time dimension of the slice
                for t in 1:size(data_slice)[3]
                    step_data = data_slice[:, :, t]
                    
                    # Apply weights & Calculate Mean
                    w_prod = step_data .* weights' # Check transpose!
                    
                    # Safe mean calculation
                    valid_mask = .!ismissing.(w_prod) .& .!isnan.(w_prod)
                    if any(valid_mask)
                        val = sum(w_prod[valid_mask]) / sum(weights'[valid_mask])
                        push!(temp_buffer, val)
                    end
                end
            end
        end
        
        # Store the year's result if we found data
        if !isempty(temp_buffer)
            push!(yearly_means, mean(temp_buffer))
            push!(valid_years, year)
            println("Year $year: $(round(mean(temp_buffer) - 273.15, digits=2))°C")
        end
    end

    # Plotting
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


visu_filtered_climatology(data_folder, weights, 1950:1960, selected_months=[1], selected_days=1:15)