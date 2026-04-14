# Chargement des fonctions :
include("function2charge.jl")
# Chargement des variables de base :
include("dataset2load.jl")

using StateSpaceModels
using CategoricalArrays
using LinearAlgebra

# Test corrélation variance/moyenne : Objectif détermination modèle additif/multiplicatif

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

#Ajouter plot

#0.05 de score de coefficiant de Pearson, on peut donc supposer l'additivité du modèle. 
# Pas de corrélation entre moyenne et variance.


using NPZ
npzwrite("month_mean_serie.npy", monthly_mean_series)
npzwrite("cycle_serie.npy", z_french_cycle)
npzwrite("ea_serie.npy", z_empirical_ea)
npzwrite("nao_serie.npy", z_empirical_nao)


# Initilisation des vecteurs (jours/mois/années)
matrix_days = compute_general_climatology(data_folder_basic, weight_prop_basic, 1950:2025, mode=:daily)
matrix_months = compute_general_climatology(data_folder_basic, weight_prop_basic, 1950:2025, mode=:monthly)
matrix_years = compute_general_climatology(data_folder_basic, weight_prop_basic, 1950:2025, mode=:yearly)

daily_mean_series = means_vector_calculation(matrix_days, weight_prop_basic)
monthly_mean_series = means_vector_calculation(matrix_months, weight_prop_basic)
annual_mean_series = means_vector_calculation(matrix_years, weight_prop_basic)

time_months = 1:Proj.length(monthly_mean_series)
time_years = 1:Proj.length(annual_mean_series)

t_norm_months = (time_months .- mean(time_months)) ./ std(time_months)
t_norm_years = (time_years .- mean(time_years)) ./ std(time_years)


df_months = DataFrame(
    Y = monthly_mean_series,
    t1 = t_norm_months,
    t2 = t_norm_months.^2
)

df_years = DataFrame(
    Y = annual_mean_series,
    t1 = t_norm_years,
    t2 = t_norm_years.^2
)





# 1. Create a repeating array of months (1 through 12) to match your data length
# Assuming your data starts in January. If it starts in another month, adjust the offset.
month_sequence = repeat(1:12, inner=1, outer=ceil(Int, length(monthly_mean_series)/12))[1:length(monthly_mean_series)]

# 2. Add it to the DataFrame and explicitly cast it as categorical
df_months.Month = categorical(month_sequence)

# 3. Rerun the regression with the seasonal dummy variables
# The GLM will automatically create 11 dummy variables to absorb the seasonal wave
poly_seasonal_model = lm(@formula(Y ~ 0 + t1 + t2 + Month), df_months)

println(poly_seasonal_model)

coefs = coef(poly_seasonal_model)

monthly_intercepts = coefs[3:14]

# The true secular intercept is the average of all 12 months
true_global_intercept = mean(monthly_intercepts)
# 3. Reconstruct ONLY the secular trend (ignoring the Month coefficients)
# coefs[1] = Intercept (January baseline)
# coefs[2] = t1 coefficient
# coefs[3] = t2 coefficient
secular_quadratic_trend = true_global_intercept .+ (coefs[1] .* df_months.t1) .+ (coefs[2] .* df_months.t2)

# 4. Detrend the raw data
# By subtracting ONLY the secular trend, the resulting array retains 
# the natural seasonality, the AMV cycle, and the weather noise.
data_for_statespace = df_months.Y .- secular_quadratic_trend
data_for_statespace= data_for_statespace .- mean(data_for_statespace)
# 5. Extract the Cycle using the highly stable State-Space architecture
final_ss_model = UnobservedComponents(
    data_for_statespace; 
    seasonal = "deterministic 12", 
    cycle = "stochastic"
)



StateSpaceModels.fit!(final_ss_model; optimizer = Optimizer(StateSpaceModels.Optim.NelderMead()), save_hyperparameter_distribution=false)
ks_final = kalman_smoother(final_ss_model)


using Plots
using LinearAlgebra



# ==============================================================================
# 1. EXTRACT THE SEASONAL HARMONIC
# ==============================================================================
# The state vector (ks_final.alpha) contains the seasonal dummies at indices 1-11.
# We multiply these states by the Observation Matrix (Z) to reconstruct the exact 
# amplitude of the seasonal wave at every time step.
# ==============================================================================
# 1. EXTRACT THE SEASONAL HARMONIC
# ==============================================================================
# In StateSpaceModels.jl, the system.Z parameter is a 1D Vector for invariant models.
Z_vec = final_ss_model.system.Z

# We isolate the 11 seasonal states (Indices 2 through 12) and take the dot product 
# against the corresponding Z weights to reconstruct the exact seasonal wave.
extracted_seasonality = [dot(Z_vec[2:12], ks_final.alpha[t][2:12]) for t in 1:length(df_months.Y)]

# ==============================================================================
# 2. EXTRACT THE NATURAL OCEANIC DRIVER
# ==============================================================================
# The cycle component begins precisely at index 13.
extracted_pure_cycle = [state[13] for state in ks_final.alpha]

extracted_trend = [state[1] for state in ks_final.alpha]

extracted_noise = data_for_statespace .- (extracted_trend .+ extracted_seasonality .+ extracted_pure_cycle)


# ==============================================================================
# 2. PLOTTING ROW 1: THE RIGID ANTHROPOGENIC BASELINE
# ==============================================================================
p1 = plot(time_months, df_months.Y, 
          label="Raw Monthly Data", 
          color=:lightgray, 
          linewidth=1,
          ylabel="Temperature (°C)",
          title="Quadratic Baseline (Accelerating Warming)",
          legend=:topleft)

# Overlay the deterministic quadratic curve from the GLM
# This uses the 'secular_quadratic_trend' array generated in the previous step
plot!(p1, time_months, secular_quadratic_trend, 
      label="Secular Trend (t + t²)", 
      color=:red, 
      linewidth=3)


# ==============================================================================
# 3. PLOTTING ROW 2: THE ISOLATED SEASONAL WAVE
# ==============================================================================
# We plot a 5-year slice (60 months) so the intra-annual shape is visible, 
# preventing the lines from compressing into a solid block of color.
zoom_slice = 1:60 
p2 = plot(time_months[zoom_slice], extracted_seasonality[zoom_slice], 
          label="12-Month Harmonic", 
          color=:forestgreen, 
          linewidth=2,
          ylabel="Amplitude",
          title="Isolated Seasonality (5-Year Zoom)")

hline!(p2, [0], color=:black, linestyle=:dash, label="")


# ==============================================================================
# 4. PLOTTING ROW 3: THE NATURAL OCEANIC DRIVER
# ==============================================================================
p3 = plot(time_months, extracted_pure_cycle, 
          label="Stochastic Cycle (C_t)", 
          color=:mediumblue, 
          linewidth=2,
          ylabel="Amplitude",
          title="Isolated Climatological Cycle (AMV/NAO)")

hline!(p3, [0], color=:black, linestyle=:dash, label="")

p4 = plot(time_months, extracted_trend, 
          label="Deterministic Offset", 
          color=:darkred, 
          linewidth=2,
          ylabel="Amplitude",
          title="Isolated Climatological Trend (Re-centering Offset)")
hline!(p4, [0], color=:black, linestyle=:dash, label="")

p5 = plot(time_months, extracted_noise, 
          label="True White Noise (ε_t)", 
          color=:black, 
          linewidth=1, # Thinner line to visualize the high-frequency scatter
          ylabel="Amplitude",
          title="Isolated Climatological Noise")
hline!(p5, [0], color=:black, linestyle=:dash, label="")

# Final Dashboard Assembly
final_dashboard = plot(p1, p2, p3, p4, p5, 
                       layout=(5, 1), 
                       size=(1000, 1200), # Increased height to accommodate 5 panels cleanly
                       margin=5Plots.mm)

display(final_dashboard)




# test DSP

using DSP



# Create deterministic seasonal dummies
month_sequence = repeat(1:12, inner=1, outer=ceil(Int, length(monthly_mean_series)/12))[1:length(monthly_mean_series)]


# Fit the Zero-Intercept GLM
poly_seasonal_model = lm(@formula(Y ~ 0 + t1 + t2 + Month), df_months)
coefs = coef(poly_seasonal_model)

# Compute True Global Baseline
monthly_intercepts = coefs[3:14]
true_global_intercept = mean(monthly_intercepts)

# Reconstruct the Secular Accelerating Trend
centered_secular_trend = true_global_intercept .+ (coefs[1] .* df_months.t1) .+ (coefs[2] .* df_months.t2)

# Reconstruct the Deterministic Seasonal Wave
# We map the 12 coefficients back to the timeline by subtracting the global mean
seasonal_coefficients_centered = monthly_intercepts .- true_global_intercept
extracted_seasonality = [seasonal_coefficients_centered[m] for m in month_sequence]

# ==============================================================================
# STEP 2: ALGEBRAIC PRE-WHITENING
# ==============================================================================
# This array is now strictly zero-mean, completely stripped of secular warming 
# and the massive 12-month seasonal variance. 
glm_residuals = df_months.Y .- (centered_secular_trend .+ extracted_seasonality)

# ==============================================================================
# STEP 3: DIGITAL SIGNAL PROCESSING (THE STOCHASTIC CYCLE)
# ==============================================================================
# Design the Butterworth Band-Pass for the AMO/NAO.
# Example: 10 years to 60 years -> 120 to 720 months.
filter_response = Bandpass(1/84, 1/24)
filter_design = Butterworth(4)

# Apply zero-phase filtering directly to the GLM residuals
extracted_oceanic_cycle = filtfilt(digitalfilter(filter_response, filter_design), glm_residuals)

# ==============================================================================
# STEP 4: ALGEBRAIC NOISE ISOLATION
# ==============================================================================
# True stochastic weather noise is whatever the DSP filter ignored
extracted_noise = glm_residuals .- extracted_oceanic_cycle


# ==============================================================================
# STEP 5: FINAL PUBLICATION DASHBOARD
# ==============================================================================
p1 = plot(time_months, df_months.Y, label="Raw Monthly Data", color=:lightgray, linewidth=1, title="Quadratic Anthropogenic Baseline")
plot!(p1, time_months, centered_secular_trend, label="Secular Trend (t + t²)", color=:red, linewidth=3)

p2 = plot(time_months[1:60], extracted_seasonality[1:60], label="12-Month Harmonic", color=:forestgreen, linewidth=2, title="Isolated Symmetrical Seasonality")
hline!(p2, [0], color=:black, linestyle=:dash, label="")

p3 = plot(time_months, extracted_oceanic_cycle, label="Stochastic Cycle (AMO)", color=:mediumblue, linewidth=2, title="DSP Extracted Multidecadal Oceanic Wave")
hline!(p3, [0], color=:black, linestyle=:dash, label="")

p4 = plot(time_months, extracted_noise, label="True White Noise (ε_t)", color=:black, linewidth=1, title="Isolated Climatological Noise")
hline!(p4, [0], color=:black, linestyle=:dash, label="")

final_dashboard = plot(p1, p2, p3, p4, layout=(4, 1), size=(1000, 1000), margin=5Plots.mm)
display(final_dashboard)







  

# Test avec données anuelles


# ==============================================================================
# 1. THE SECULAR ACCELERATION (GLM)
# ==============================================================================
# We revert to a standard intercept because seasonal masking no longer exists.
poly_model = lm(@formula(Y ~ 1 + t1 + t2), df_years)
println(poly_model)

coefs = coef(poly_model)

# coefs[1] = Intercept (True Annual Global Baseline)
# coefs[2] = t1 coefficient
# coefs[3] = t2 coefficient
secular_quadratic_trend = coefs[1] .+ (coefs[2] .* df_years.t1) .+ (coefs[3] .* df_years.t2)

# Detrend the raw data
data_for_statespace = df_years.Y .- secular_quadratic_trend

# ==============================================================================
# 2. STATE-SPACE EXTRACTION
# ==============================================================================
# The seasonal parameter is strictly omitted.
final_ss_model = UnobservedComponents(
    data_for_statespace; 
    cycle = "stochastic"
)

StateSpaceModels.fit!(
    final_ss_model; 
    optimizer = Optimizer(StateSpaceModels.Optim.NelderMead()), 
    save_hyperparameter_distribution=false
)

ks_final = kalman_smoother(final_ss_model)

# ==============================================================================
# 3. STATE VECTOR SLICING & NOISE ISOLATION
# ==============================================================================
# The state vector now has exactly 3 elements.
# Index 1 = Local Level (residual mean offset)
# Index 2 = Stochastic Cycle (C_t)
# Index 3 = Conjugate Cycle (C_t*)

extracted_trend_offset = [state[1] for state in ks_final.alpha]
extracted_pure_cycle   = [state[2] for state in ks_final.alpha]

# Calculate True White Noise (ε_t) via algebraic subtraction
extracted_noise = data_for_statespace .- (extracted_trend_offset .+ extracted_pure_cycle)

# ==============================================================================
# 4. DASHBOARD RENDERING (3-PANEL LAYOUT)
# ==============================================================================

p1 = plot(time_years, df.Y, 
          label="Raw Annual Data", 
          color=:lightgray, 
          linewidth=2,
          ylabel="Temperature (°C)",
          title="Quadratic Baseline (Accelerating Warming)",
          legend=:topleft)

# Overlay the deterministic quadratic curve
plot!(p1, time_years, secular_quadratic_trend, 
      label="Secular Trend (t + t²)", 
      color=:red, 
      linewidth=3)

p2 = plot(time_years, extracted_pure_cycle, 
          label="Stochastic Cycle (C_t)", 
          color=:mediumblue, 
          linewidth=2,
          ylabel="Amplitude",
          title="Isolated Climatological Cycle (AMV/NAO)")
hline!(p2, [0], color=:black, linestyle=:dash, label="")

p3 = plot(time_years, extracted_noise, 
          label="True White Noise (ε_t)", 
          color=:black, 
          linewidth=1,
          ylabel="Amplitude",
          title="Isolated Inter-annual Noise")
hline!(p3, [0], color=:black, linestyle=:dash, label="")

# Final Dashboard Assembly
# We condense to a 3-panel layout, omitting the seasonal and dead-trend panels
final_dashboard = plot(p1, p2, p3, 
                       layout=(3, 1), 
                       size=(1000, 900), 
                       margin=5Plots.mm)

display(final_dashboard)









using StatsBase
using CSV
NAO_df = CSV.read("demo_france/NAO.txt", DataFrame, ignorerepeated=true, delim=' ', header=false)
NAO_serie = NAO_df.Column3[1:912]

plot(time_months, NAO_serie)

# ==============================================================================
# 2. ENFORCE SPECTRAL CONSISTENCY
# ==============================================================================
# You MUST use the exact same filter bounds here that you used on the French data.
# The example below assumes you targeted the inter-annual NAO (24 to 84 months).
# If you targeted the decadal NAO (84 to 120 months), adjust the Bandpass accordingly.
filter_response = Bandpass(1/84, 1/24)
filter_design = Butterworth(4)

# Apply zero-phase filtering to the empirical NAO index
filtered_empirical_nao = filtfilt(digitalfilter(filter_response, filter_design), NAO_serie)


function standardize(array)
    return (array .- mean(array)) ./ std(array)
end

z_french_cycle = standardize(extracted_oceanic_cycle)
z_empirical_nao = standardize(filtered_empirical_nao)

# ==============================================================================
# 4. CROSS-CORRELATION FUNCTION (CCF)
# ==============================================================================
# Define the maximum lag to test. 60 months (5 years) is standard for inter-annual cycles.
max_lag = 300 

# Compute Pearson correlation across phase shifts.
# Order matters: crosscor(x, y) tests if 'x' leads 'y' at positive lags.
ccf_values = crosscor(z_empirical_nao, z_french_cycle, -max_lag:max_lag)

# Define strict 95% statistical significance boundaries
n_samples = length(z_french_cycle)
ci_bound = 1.96 / sqrt(n_samples)

# ==============================================================================
# 5. VISUAL VALIDATION DASHBOARD
# ==============================================================================
time_axis = 1:n_samples
lag_axis = -max_lag:max_lag

# Panel 1: Phase Synchronization Overlay
p1 = plot(time_axis, z_empirical_nao, 
          label="Filtered Empirical NAO Index", 
          color=:gray, 
          linewidth=2,
          ylabel="Z-Score",
          title="Spectral Synchronization Overlay (Matched Frequencies)")

plot!(p1, time_axis, z_french_cycle, 
      label="Extracted French Cycle", 
      color=:blue, 
      linewidth=2)

hline!(p1, [0], color=:black, linestyle=:dash, label="")

# Panel 2: The Correlogram
p2 = bar(lag_axis, ccf_values, 
         label="Correlation Coefficient (r)", 
         color=:darkblue, 
         linecolor=:transparent,
         xlabel="Lag (Months)", 
         ylabel="Correlation (r)",
         title="Cross-Correlation Function")

# Add significance boundaries
hline!(p2, [ci_bound], color=:red, linestyle=:dash, label="95% Significance Limit")
hline!(p2, [-ci_bound], color=:red, linestyle=:dash, label="")
hline!(p2, [0], color=:black, linewidth=1, label="")

validation_dashboard = plot(p1, p2, layout=(2, 1), size=(1000, 800), margin=5Plots.mm)
display(validation_dashboard)


# Big correlation at 135 months (11-12 years).
# Testons une modélisation.




using RData # Optional, depending on how you score R²

# Assume 'z_french_cycle' and 'z_empirical_nao' are your arrays from the previous step.

# ==============================================================================
# 1. BUILD THE DUAL-FORCING DATAFRAME
# ==============================================================================
df_dual = DataFrame(
    French_Temp = z_french_cycle,
    NAO_Lag0    = z_empirical_nao
)

# ==============================================================================
# 2. GENERATE THE OCEANIC MEMORY VECTOR
# ==============================================================================
lag_ocean = 80

# Shift the NAO array forward by 135 months, padding the start with missing values
df_dual.NAO_Lag135 = [fill(missing, lag_ocean); z_empirical_nao[1:end-lag_ocean]]

# Drop the first 11.25 years to align the temporal matrices
df_model = dropmissing(df_dual)

# ==============================================================================
# 3. FIT THE DISTRIBUTED LAG MODEL
# ==============================================================================
# We use the zero-intercept architecture (0 + ...)
dual_forcing_model = lm(@formula(French_Temp ~ 0 + NAO_Lag0 + NAO_Lag135), df_model)

# ==============================================================================
# 4. DIAGNOSTIC OUTPUT
# ==============================================================================
println(dual_forcing_model)
println("Total Variance Explained (R²): ", r2(dual_forcing_model))

# Generate the predictive model's internal timeline
predicted_french_climate = predict(dual_forcing_model)

time_axis_model = 81:912

# Extract the actual values the model was trying to predict directly from the clean dataframe
actual_french_climate = df_model.French_Temp

# ==============================================================================
# 2. RENDER THE PREDICTIVE VALIDATION DASHBOARD
# ==============================================================================
# Instantiate the plot and assign it to a variable (p_dlm)
p_dlm = plot(time_axis_model, actual_french_climate, 
          label="Actual Extracted French Cycle (C_t)", 
          color=:gray, 
          linewidth=2,
          ylabel="Z-Score",
          title="DLM Validation: Dual-Forcing Oceanic Prediction")

# Overlay the predictive array generated by the 0 + 135 Lag Model
plot!(p_dlm, time_axis_model, predicted_french_climate, 
      label="DLM Prediction (Kinematic + Thermodynamic)", 
      color=:red, 
      linewidth=2)

# Add the strict zero-mean baseline
hline!(p_dlm, [0], color=:black, linestyle=:dash, label="")

# Display the final output
display(p_dlm)





# Test avec EA (east atlantic pattern)

EA_df = CSV.read("demo_france/EA.txt", DataFrame, ignorerepeated=true, delim=' ', header=1)
EA_serie = EA_df.INDEX[1:912]

plot(time_months, EA_serie)

# ==============================================================================
# 2. ENFORCE SPECTRAL CONSISTENCY
# ==============================================================================
# You MUST use the exact same filter bounds here that you used on the French data.
# The example below assumes you targeted the inter-annual NAO (24 to 84 months).
# If you targeted the decadal NAO (84 to 120 months), adjust the Bandpass accordingly.
filter_response = Bandpass(1/84, 1/24)
filter_design = Butterworth(4)

# Apply zero-phase filtering to the empirical NAO index
filtered_empirical_ea = filtfilt(digitalfilter(filter_response, filter_design), EA_serie)


function standardize(array)
    return (array .- mean(array)) ./ std(array)
end

z_french_cycle = standardize(extracted_oceanic_cycle)
z_empirical_ea = standardize(filtered_empirical_ea)

# ==============================================================================
# 4. CROSS-CORRELATION FUNCTION (CCF)
# ==============================================================================
# Define the maximum lag to test. 60 months (5 years) is standard for inter-annual cycles.
max_lag = 300 

# Compute Pearson correlation across phase shifts.
# Order matters: crosscor(x, y) tests if 'x' leads 'y' at positive lags.
ccf_values = crosscor(z_empirical_ea, z_french_cycle, -max_lag:max_lag)

# Define strict 95% statistical significance boundaries
n_samples = length(z_french_cycle)
ci_bound = 1.96 / sqrt(n_samples)

# ==============================================================================
# 5. VISUAL VALIDATION DASHBOARD
# ==============================================================================
time_axis = 1:n_samples
lag_axis = -max_lag:max_lag

# Panel 1: Phase Synchronization Overlay
p1 = plot(time_axis, z_empirical_ea, 
          label="Filtered Empirical EA Index", 
          color=:gray, 
          linewidth=2,
          ylabel="Z-Score",
          title="Spectral Synchronization Overlay (Matched Frequencies)")

plot!(p1, time_axis, z_french_cycle, 
      label="Extracted French Cycle", 
      color=:blue, 
      linewidth=2)

hline!(p1, [0], color=:black, linestyle=:dash, label="")

# Panel 2: The Correlogram
p2 = bar(lag_axis, ccf_values, 
         label="Correlation Coefficient (r)", 
         color=:darkblue, 
         linecolor=:transparent,
         xlabel="Lag (Months)", 
         ylabel="Correlation (r)",
         title="Cross-Correlation Function")

# Add significance boundaries
hline!(p2, [ci_bound], color=:red, linestyle=:dash, label="95% Significance Limit")
hline!(p2, [-ci_bound], color=:red, linestyle=:dash, label="")
hline!(p2, [0], color=:black, linewidth=1, label="")

validation_dashboard = plot(p1, p2, layout=(2, 1), size=(1000, 800), margin=5Plots.mm)
display(validation_dashboard)



df_model.EA_Lag0 = z_empirical_ea[136:end] # Adjust indices strictly to match df_model length

# The Master Multivariate DLM
master_climate_model = lm(@formula(French_Temp ~ 0 + NAO_Lag0 + NAO_Lag135 + EA_Lag0), df_model)
master_prediction = predict(master_climate_model)

p_dlm = plot(time_axis_model, actual_french_climate, 
          label="Actual Extracted French Cycle (C_t)", 
          color=:gray, 
          linewidth=2,
          ylabel="Z-Score",
          title="DLM Validation: Dual-Forcing Oceanic Prediction")

# Overlay the predictive array generated by the 0 + 135 Lag Model
plot!(p_dlm, time_axis_model, master_prediction, 
      label="DLM Prediction (Kinematic + Thermodynamic)", 
      color=:red,
      linewidth=2)

# Add the strict zero-mean baseline
hline!(p_dlm, [0], color=:black, linestyle=:dash, label="")

# Display the final output
display(p_dlm)



compound_model = lm(
    @formula(French_Temp ~ 0 + NAO_Lag0 + NAO_Lag135 + EA_Lag0 + NAO_Lag0&EA_Lag0 + NAO_Lag135&EA_Lag0), 
    df_model
)

compound_prediction = predict(compound_model)

p_dlm = plot(time_axis_model, actual_french_climate, 
          label="Actual Extracted French Cycle (C_t)", 
          color=:gray, 
          linewidth=2,
          ylabel="Z-Score",
          title="DLM Validation: Dual-Forcing Oceanic Prediction")

# Overlay the predictive array generated by the 0 + 135 Lag Model
plot!(p_dlm, time_axis_model, compound_prediction, 
      label="DLM Prediction (Kinematic + Thermodynamic)", 
      color=:red,
      linewidth=2)

# Add the strict zero-mean baseline
hline!(p_dlm, [0], color=:black, linestyle=:dash, label="")

# Display the final output
display(p_dlm)

# Echec avec relations non linéaires. Meilleur modèle en additif.





df_master = DataFrame(
    French_Temp = z_french_cycle,
    NAO_Lag0    = z_empirical_nao,
    EA_Lag0     = z_empirical_ea
)


lag_ocean = 135
df_master.NAO_Lag135 = [fill(missing, lag_ocean); z_empirical_nao[1:end-lag_ocean]]

# Shift NAO by 1 month to establish the previous state
df_master.NAO_Lag1 = [missing; z_empirical_nao[1:end-1]]

df_master.NAO_Momentum = df_master.NAO_Lag0 .- df_master.NAO_Lag1

df_ortho = dropmissing(df_master)




# ==============================================================================
# 2. FIT THE ORTHOGONALIZED MODEL
# ==============================================================================
# We inject NAO_Momentum instead of NAO_Lag1 to prevent the matrix collapse
momentum_model = lm(
    @formula(French_Temp ~ 0 + NAO_Lag0 + NAO_Momentum + NAO_Lag135 + EA_Lag0), 
    df_ortho
)

momentum_prediction = predict(momentum_model)

p_dlm = plot(time_axis_model, actual_french_climate, 
          label="Actual Extracted French Cycle (C_t)", 
          color=:gray, 
          linewidth=2,
          ylabel="Z-Score",
          title="DLM Validation: Dual-Forcing Oceanic Prediction")

# Overlay the predictive array generated by the 0 + 135 Lag Model
plot!(p_dlm, time_axis_model, momentum_prediction, 
      label="DLM Prediction (Kinematic + Thermodynamic)", 
      color=:red,
      linewidth=2)

# Add the strict zero-mean baseline
hline!(p_dlm, [0], color=:black, linestyle=:dash, label="")

# Display the final output
display(p_dlm)


# Effet positif clair du modèle.


# Create asymmetric vectors for the EA pattern
df_ortho.EA_Positive = max.(df_ortho.EA_Lag0, 0.0)
df_ortho.EA_Negative = min.(df_ortho.EA_Lag0, 0.0)

# Run the asymmetric model
asymmetric_model = lm(
    @formula(French_Temp ~ 0 + NAO_Lag0 + NAO_Momentum + NAO_Lag135 + EA_Positive + EA_Negative), 
    df_ortho
)


# EA_Positive n'a pas d'effet sur la prédiction. On peut donc l'enlever.


master_asymmetric_model = lm(
    @formula(French_Temp ~ 0 + NAO_Lag0 + NAO_Momentum + NAO_Lag135 + EA_Negative), 
    df_ortho
)


master_asymmetric_prediction = predict(master_asymmetric_model)

p_dlm = plot(time_axis_model, actual_french_climate, 
          label="Actual Extracted French Cycle (C_t)", 
          color=:gray, 
          linewidth=2,
          ylabel="Z-Score",
          title="DLM Validation: Dual-Forcing Oceanic Prediction")

# Overlay the predictive array generated by the 0 + 135 Lag Model
plot!(p_dlm, time_axis_model, master_asymmetric_prediction, 
      label="DLM Prediction (Kinematic + Thermodynamic)", 
      color=:red,
      linewidth=2)

plot!(p_dlm, time_axis_model, momentum_prediction, 
      label="DLM Prediction (Kinematic + Thermodynamic)", 
      color=:green,
      linewidth=2)


# Add the strict zero-mean baseline
hline!(p_dlm, [0], color=:black, linestyle=:dash, label="")

# Display the final output
display(p_dlm)




# Testons l'influence des différentes composantes au cours du temps.



window_size = 360 # 30 years * 12 months
n_months = nrow(df_ortho)
ea_coefficients = Float64[]
a135_coefficents = Float64[]
a1_coefficients = Float64[]
diff_coefficients = Float64[]

# Slide a 30-year window across the timeline
for start_idx in 1:(n_months - window_size)
    end_idx = start_idx + window_size - 1
    window_df = df_ortho[start_idx:end_idx, :]
    
    window_model = lm(
        @formula(French_Temp ~ 0 + NAO_Lag0 + NAO_Momentum + NAO_Lag135 + EA_Lag0), 
        window_df
    )
    
    # Extract just the EA coefficient for this 30-year era
    push!(ea_coefficients, coef(window_model)[4]) 
    push!(a135_coefficents, coef(window_model)[3]) 
    push!(a1_coefficients, coef(window_model)[1]) 
    push!(diff_coefficients, coef(window_model)[2]) 
end

plot(1:417, ea_coefficients)
plot!(1:417, a135_coefficents)
plot!(1:417, a1_coefficients)
plot!(1:417, diff_coefficients)

# On peut clairement voir que EA et Diff sont changeants au cours du temps, 
# ce qui n'est pas le cas de a1 et a135.

# Testons en séparant la série au mois 

# ==============================================================================
# 1. DEFINE THE STRUCTURAL BREAK DUMMIES
# ==============================================================================
# Assuming your df_ortho has the :Year column from the initial CSV parse.
# If not, you can use row indices: (1:nrow(df_ortho)) .< 350
historical_dummy = (1:nrow(df_ortho)) .< 350
modern_dummy = (1:nrow(df_ortho)) .>= 350

# ==============================================================================
# 2. CREATE THE REGIME-SPECIFIC VECTORS
# ==============================================================================
# We split the Momentum vector
df_ortho.Momentum_Historical = df_ortho.NAO_Momentum .* historical_dummy
df_ortho.Momentum_Modern     = df_ortho.NAO_Momentum .* modern_dummy

# We split the East Atlantic vector
df_ortho.EA_Historical = df_ortho.EA_Lag0 .* historical_dummy
df_ortho.EA_Modern     = df_ortho.EA_Lag0 .* modern_dummy

# ==============================================================================
# 3. FIT THE STRUCTURAL REGIME MODEL
# ==============================================================================
regime_model = lm(
    @formula(French_Temp ~ 0 + NAO_Lag0 + NAO_Lag135 + 
             Momentum_Historical + Momentum_Modern + 
             EA_Historical + EA_Modern), 
    df_ortho
)
regime_prediction = predict(regime_model)

p_dlm = plot(time_axis_model, actual_french_climate, 
          label="Actual Extracted French Cycle (C_t)", 
          color=:gray, 
          linewidth=2,
          ylabel="Z-Score",
          title="DLM Validation: Dual-Forcing Oceanic Prediction")

# Overlay the predictive array generated by the 0 + 135 Lag Model
plot!(p_dlm, time_axis_model, master_asymmetric_prediction, 
      label="DLM Prediction (Kinematic + Thermodynamic)", 
      color=:red,
      linewidth=2)

plot!(p_dlm, time_axis_model, momentum_prediction, 
      label="DLM Prediction (Kinematic + Thermodynamic)", 
      color=:green,
      linewidth=2)

plot!(p_dlm, time_axis_model, regime_prediction, 
      label="DLM Prediction (Kinematic + Thermodynamic)", 
      color=:blue,
      linewidth=2)




# Test avec les variations saisonnières (été et hiver)
# A refaire plus précisement.

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