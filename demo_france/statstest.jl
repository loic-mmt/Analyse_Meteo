# Chargement des fonctions :
include("function2charge.jl")
# Chargement des variables de base :
include("dataset2load.jl")


# Determination 

vect_vars = []
vect_means = []
for index in 1950:2025
    brut_matrix = compute_general_climatology(data_folder_basic, weight_prop_basic, index, mode=:hourly)
    vector = means_vector_calculation(brut_matrix, weight_prop_basic)
    mean_temp = mean(vector)
    variance = var(vector)
    push!(vect_vars, variance)
    push!(vect_means, mean_temp)
end
cor(vect_vars, vect_means)

plot()

#0.05 de score de coefficiant de Pearson, on peut donc supposer l'additivité du modèle. 
# Pas de corrélation entre moyenne et variance.



matrix_days = compute_general_climatology(data_folder_basic, weight_prop_basic, 1950:2025, mode=:daily)
matrix_months = compute_general_climatology(data_folder_basic, weight_prop_basic, 1950:2025, mode=:monthly)
daily_mean_series = means_vector_calculation(matrix_days, weight_prop_basic)
monthly_mean_series = means_vector_calculation(matrix_months, weight_prop_basic)


using StateSpaceModels
using SeasonalTrendLoess

res_monthly_stl = stl(monthly_mean_series, 12; robust=true, ns=999)

# Create the strictly non-seasonal dataset
# This array now contains ONLY the long-term trend, the cycles, and the noise
true_deseasonalized_data = monthly_mean_series - res_monthly_stl.seasonal

# STEP 2: The Low-Dimensional State-Space Model
# CRITICAL: Omit the 'seasonal' parameter entirely. 
hybrid_model = UnobservedComponents(
    true_deseasonalized_data, 
    trend="local linear trend", 
    cycle="stochastic"
)

StateSpaceModels.fit!(hybrid_model; optimizer = Optimizer(StateSpaceModels.Optim.NelderMead()))
ks = kalman_smoother(hybrid_model)



trend_level = [state[1] for state in ks.alpha]
trend_slope = [state[2] for state in ks.alpha]
stochastic_cycle = [state[3] for state in ks.alpha]
trend_variance = [cov_matrix[1, 1] for cov_matrix in ks.V]
trend_std = sqrt.(trend_variance)
upper_bound = trend_level .+ (1.96 .* trend_std)
lower_bound = trend_level .- (1.96 .* trend_std)

time_index = 1:length(trend_level)

# ---------------------------------------------------------
# Panel 1: The Global Trend and Uncertainty Overlaid on Data
# ---------------------------------------------------------
# First, plot the raw or deseasonalized data as a faint background layer
p1 = plot(time_index, monthly_mean_series, 
          label="Deseasonalized Data", 
          color=:lightgray, 
          linewidth=1,
          ylabel="Value",
          title="Global Climatic Trend (μ_t)")

# Next, overlay the smoothed trend line and use the `ribbon` argument
# to automatically generate the 95% confidence interval shaded bands
plot!(p1, time_index, trend_level, 
      ribbon=(1.96 .* trend_std, 1.96 .* trend_std), 
      fillalpha=0.3, # Sets the transparency of the confidence band
      label="Smoothed Trend ± 95% CI", 
      color=:blue, 
      linewidth=2)


# ---------------------------------------------------------
# Panel 2: The Instantaneous Slope (Climatic Drift)
# ---------------------------------------------------------
p2 = plot(time_index, trend_slope,
          label="Instantaneous Slope (β_t)",
          color=:darkorange,
          linewidth=2,
          ylabel="Rate of Change",
          title="Climatic Drift")

# Add a horizontal reference line at zero. 
# If the line is above zero, the climatic trend is accelerating upwards.
hline!(p2, [0], color=:black, linestyle=:dash, label="")


# ---------------------------------------------------------
# Panel 3: The Isolated Stochastic Cycle
# ---------------------------------------------------------
p3 = plot(time_index, stochastic_cycle, 
          label="Stochastic Cycle (C_t)", 
          color=:forestgreen, 
          linewidth=2,
          ylabel="Amplitude",
          title="Glocal Cyclical Phenomena")

# Add a horizontal zero line to clearly identify cyclical peaks and troughs
hline!(p3, [0], color=:black, linestyle=:dash, label="")


# ---------------------------------------------------------
# Final Assembly
# ---------------------------------------------------------
# Concatenate the three panels vertically
final_dashboard = plot(p1, p2, p3, 
                       layout=(3, 1), 
                       size=(900, 900), 
                       link=:x, # Crucial: links the x-axes so zooming scales all panels simultaneously
                       left_margin=5Plots.mm)

# Render the plot to the display
display(final_dashboard)


# Test avec les variations saisonnières (été et hiver)

summer_matrix = compute_general_climatology(data_folder_basic, weight_prop_basic, 1950:2025, selected_months= [7,8], mode=:yearly)
summer_means = means_vector_calculation(summer_matrix, weight_prop_basic)
summer_rolling_var = runvar(summer_means, 30)[30:end]



winter_matrix = compute_general_climatology(data_folder_basic, weight_prop_basic, 1950:2025, selected_months= [1,2], mode=:yearly)
winter_means = means_vector_calculation(winter_matrix, weight_prop_basic)
winter_rolling_var = runvar(winter_means, 30)[30:end]


plot(1979:2025, summer_rolling_var, 
    label="summer variance evolution",
    color=:orange, 
    )

plot!(1979:2025, winter_rolling_var, 
    label="winter variance evolution",
    color=:blue, 
    )


# modélisation pour été et hiver.

winter_model = UnobservedComponents(winter_means, trend="local linear trend", cycle="stochastic")
StateSpaceModels.fit!(winter_model; optimizer = Optimizer(StateSpaceModels.Optim.NelderMead()))
ks_w = kalman_smoother(winter_model)

# Extract the Winter Cycle (Index 3: 1=Level, 2=Slope, 3=Cycle)
winter_cycle = [state[3] for state in ks_w.alpha]

# ---------------------------------------------------------
# 2. Independent Extraction for Summer
# ---------------------------------------------------------
summer_model = UnobservedComponents(summer_means, trend="local linear trend", cycle="stochastic")
StateSpaceModels.fit!(summer_model; optimizer = Optimizer(StateSpaceModels.Optim.NelderMead()))
ks_s = kalman_smoother(summer_model)

# Extract the Summer Cycle
summer_cycle = [state[3] for state in ks_s.alpha]

# ---------------------------------------------------------
# 3. Mathematical Correlation
# ---------------------------------------------------------
# Calculate the Pearson correlation coefficient between the two cycles
cycle_correlation = cor(winter_cycle, summer_cycle)
println("Pearson Correlation between Winter and Summer cycles: ", cycle_correlation)

# ---------------------------------------------------------
# 4. Visual Comparison
# ---------------------------------------------------------
# Plot both cycles on the same axis to visually inspect phase alignment and amplitude divergence
plot(1950:2025, winter_cycle, 
     label="Winter Cycle (C_t)", 
     color=:blue, 
     linewidth=2,
     title="Seasonal Cyclical Divergence",
     ylabel="Cycle Amplitude")

plot!(1950:2025, summer_cycle, 
      label="Summer Cycle (C_t)", 
      color=:orange, 
      linewidth=2)

hline!([0], color=:black, linestyle=:dash, label="")

years = 1950:2025

# ==============================================================================
# 1. DATA EXTRACTION (WINTER)
# ==============================================================================
w_level = [state[1] for state in ks_w.alpha]
w_slope = [state[2] for state in ks_w.alpha]
w_cycle = [state[3] for state in ks_w.alpha]

# Extract variance for the Winter trend to calculate 95% Confidence Intervals
w_trend_var = [cov_matrix[1, 1] for cov_matrix in ks_w.V]
w_std = sqrt.(w_trend_var)

# ==============================================================================
# 2. DATA EXTRACTION (SUMMER)
# ==============================================================================
s_level = [state[1] for state in ks_s.alpha]
s_slope = [state[2] for state in ks_s.alpha]
s_cycle = [state[3] for state in ks_s.alpha]

# Extract variance for the Summer trend to calculate 95% Confidence Intervals
s_trend_var = [cov_matrix[1, 1] for cov_matrix in ks_s.V]
s_std = sqrt.(s_trend_var)

# ==============================================================================
# 3. PLOTTING ROW 1: THE GLOBAL TRENDS (μ_t)
# ==============================================================================
# Winter Trend Panel
p1 = plot(years, winter_means, label="Raw Winter Mean", color=:lightgray, linewidth=1, title="Winter Trend (μ_t)", ylabel="Temperature (°C)", legend=:topleft)
plot!(p1, years, w_level, ribbon=(1.96 .* w_std, 1.96 .* w_std), fillalpha=0.3, label="Trend ± 95% CI", color=:blue, linewidth=2)

# Summer Trend Panel
p2 = plot(years, summer_means, label="Raw Summer Mean", color=:lightgray, linewidth=1, title="Summer Trend (μ_t)", legend=:topleft)
plot!(p2, years, s_level, ribbon=(1.96 .* s_std, 1.96 .* s_std), fillalpha=0.3, label="Trend ± 95% CI", color=:orange, linewidth=2)


# ==============================================================================
# 4. PLOTTING ROW 2: THE CLIMATIC DRIFT / SLOPE (β_t)
# ==============================================================================
# Winter Slope Panel
p3 = plot(years, w_slope, label="Winter Slope (β_t)", color=:darkblue, linewidth=2, title="Winter Drift", ylabel="Rate of Change")
hline!(p3, [0], color=:black, linestyle=:dash, label="")

# Summer Slope Panel
p4 = plot(years, s_slope, label="Summer Slope (β_t)", color=:darkorange, linewidth=2, title="Summer Drift")
hline!(p4, [0], color=:black, linestyle=:dash, label="")


# ==============================================================================
# 5. PLOTTING ROW 3: THE NATURAL CYCLE (C_t)
# ==============================================================================
# Winter Cycle Panel
p5 = plot(years, w_cycle, label="Winter Cycle (AMV)", color=:cyan, linewidth=2, title="Winter Volatility", ylabel="Cycle Amplitude")
hline!(p5, [0], color=:black, linestyle=:dash, label="")

# Summer Cycle Panel
p6 = plot(years, s_cycle, label="Summer Cycle (AMV)", color=:gold, linewidth=2, title="Summer Volatility")
hline!(p6, [0], color=:black, linestyle=:dash, label="")


# ==============================================================================
# 6. FINAL DASHBOARD ASSEMBLY
# ==============================================================================
# Define a 3x2 layout grid and map the individual plots into it
dashboard = plot(p1, p2, 
                 p3, p4, 
                 p5, p6, 
                 layout=(3, 2), 
                 size=(1200, 1000), 
                 link=:x, # Links the x-axes for synchronized zooming
                 margin=5Plots.mm)

display(dashboard)

total_rolling_var = runvar(stochastic_cycle, 360) # 360 months = 30 years

# 2. Extract only the valid saturated data (matching your 1979-2025 window)
# We must aggregate the monthly total_rolling_var to yearly to match the winter array
# Assuming you extract the July value of total_rolling_var to represent the year's center
valid_total_var = total_rolling_var[12*29+7 : 12 : end] 

# 3. Construct a DataFrame with the two valid variance arrays
df_variance = DataFrame(Total_Cycle_Var = valid_total_var, Winter_Var = winter_rolling_var)

# 4. Execute the OLS Regression
variance_model = lm(@formula(Total_Cycle_Var ~ Winter_Var), df_variance)

# 5. Output the results
println(variance_model)
println("R-squared: ", r2(variance_model))



# Verification du cycle avec les données officielles (NAO)



using StatsBase

# Assuming you have loaded a standard climate index (e.g., ONI) for the exact same 900 months
# historical_enso_index = [...] 

# Calculate the cross-correlation for lags between -24 and +24 months
max_lag = 24
lags = -max_lag:max_lag

# Compute the cross-correlation array
# Note: stochastic_cycle is the C_t array extracted from your Kalman Smoother
ccf_values = crosscor(historical_enso_index, stochastic_cycle, lags)

# Find the lag with the highest absolute correlation
optimal_lag_index = argmax(abs.(ccf_values))
optimal_lag_months = lags[optimal_lag_index]
peak_correlation = ccf_values[optimal_lag_index]

println("Maximum correlation of $peak_correlation found at a lag of $optimal_lag_months months.")