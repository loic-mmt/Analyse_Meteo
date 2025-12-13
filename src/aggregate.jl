using NCDatasets
using Glob
using Statistics
using Plots

data_folder = "Analyse_Meteo/data/raw-yearly-combined/era5_fr_t2m"
weight_file = "Analyse_Meteo/loic/weights_france_final.nc"
variable_name = "t2m"


println("Loading weights...")
ds_w = NCDataset(weight_file)
# Load weights into memory (2D array: lon x lat)
# Replace 'weights' with the actual variable name in your mask file
weights = ds_w["final_weights"][:,:] 
close(ds_w)

# Calculate the sum of weights (Denominator for weighted average)
# We assume 0.0 in the weights means "outside the region"
total_weight = sum(weights)

# ==========================================
# 3. Iterate and Calculate
# ==========================================
# Find all .nc files (e.g., era5_2010_01.nc, era5_2010_02.nc...)
files = glob("*.nc", data_folder)


println("Processing $(length(files)) files...")


function calculate_climatology(files, weights, variable_name)
    global_sum = 0.0
    count_steps = 0
    for file in files
    # 1. Open the dataset manually
        ds = NCDataset(file)
    
    # --- Start of Processing ---
    # Get handle to the variable
        var = ds[variable_name]
    
        n_time = size(var)[3]
        for t in 1:n_time
        # Read ONE time slice
            data_slice = var[:,:,t]
        
            weighted_product = data_slice .* weights'
            valid_mask = .!ismissing.(weighted_product) .& .!isnan.(weighted_product)
            current_sum = sum(weighted_product[valid_mask])
            current_total_weight = sum(weights'[valid_mask])

            # 5. Calculate true weighted mean
            spatial_mean = current_sum / current_total_weight
        
        # Accumulate
            global_sum += spatial_mean
            count_steps += 1
        end
    # --- End of Processing ---

    # 2. Close the dataset manually (Important!)
        close(ds)
    
        println("Finished file: $file")
    end
    return(global_sum/count_steps)
end
years = 1950:1960
temp_mean = Float64[]
for year in years
    println("Processing Year: $year")
    
    # 1. Select files for just this year
    files_year = glob("*$(year)*.nc", data_folder)
    
    # 2. Run your calculation function on just these files
    # (Using the function we wrote in the previous step)
    mean_value = calculate_climatology(files_year, weights, total_weight, variable_name)
    push!(temp_mean, mean_value)
    #mean <- mean_value/step
    println("Mean for $year: $mean_value")
end
temp_mean_celcius = temp_mean .- 273.15
plot(years, temp_mean_celcius)


function visu_temp_mean(data_folder::String, weights::Matrix{Float64}, year_range, month_range=(1,2,3,4,5,6,7,8,9,10,11,12), variable_name="t2m")
    temp_mean = Float64[]
    for year in year_range
        println("Processing Year: $year")
        temp_mean_month = Float64[]
        for month in month_range
            month = lpad(month, 2, '0')
            println("Processing Month: $month")

            # 1. Select files for just this month
            files = glob("*$year_$month*.nc", data_folder)
    
            # 2. Run your calculation function on just these files
            # (Using the function we wrote in the previous step)
            mean_value_month = calculate_climatology(files, weights, variable_name)
            push!(temp_mean_month, mean_value_month)
        end
        mean_value = mean(temp_mean_month)
        push!(temp_mean, mean_value)
        println("Mean for $year: $mean_value")
    end
    temp_mean_celcius = temp_mean .- 273.15
    plot(time, temp_mean_celcius)
    return(temp_mean_celcius)
end

visu_temp_mean(data_folder, weights, years)