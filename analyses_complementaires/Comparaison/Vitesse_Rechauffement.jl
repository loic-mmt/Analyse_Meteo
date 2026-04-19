# Vitesse de rechauffement ET COMPARAISON (FRANCE VS CANADA)

using Statistics
using Distributions
using Printf
using Plots

# 1. CONFIGURATION ET CHEMINS
project_dir = abspath(joinpath(@__DIR__, "..", ".."))
include(joinpath(project_dir, "demo_france", "function2charge.jl"))
include(joinpath(project_dir, "demo_france", "dataset2load.jl"))

cd(project_dir)

plot_dir = joinpath(@__DIR__, "plot")
isdir(plot_dir) || mkpath(plot_dir)

# 2. FONCTIONS MATHÉMATIQUES 

# Équivalent de scipy.stats.linregress
function linear_trend(x, y)
    n = length(x)
    x_mean = mean(x)
    y_mean = mean(y)
    
    ss_xx = sum((x .- x_mean).^2)
    ss_xy = sum((x .- x_mean) .* (y .- y_mean))
    
    slope = ss_xy / ss_xx
    intercept = y_mean - slope * x_mean
    pred = intercept .+ slope .* x
    
    ss_tot = sum((y .- y_mean).^2)
    ss_res = sum((y .- pred).^2)
    r2 = 1 - (ss_res / ss_tot)
    
    df = n - 2
    mse = ss_res / df
    se = sqrt(mse / ss_xx)
    
    t_stat = slope / se
    p_val = 2 * (1 - cdf(TDist(df), abs(t_stat)))
    
    return slope, intercept, r2, p_val, se, pred
end

# Équivalent de scipy.ndimage.uniform_filter1d (Lissage)
function rolling_mean(x, w)
    half = div(w, 2)
    n = length(x)
    return [mean(x[max(1, i-half):min(n, i+half)]) for i in 1:n]
end

# 3. CHARGEMENT ET CALCUL DES ANOMALIES (Réf: 1960-1990)

println("Chargement des données et calcul des moyennes spatiales...")

YEARS = 1950:2025
ref_start, ref_end = 1960, 1990
idx_ref = findall(y -> y >= ref_start && y <= ref_end, YEARS)

# -- FRANCE --
m_fr = compute_general_climatology(data_folder_basic, weight_prop_basic, YEARS, mode=:yearly)
temp_fr_abs = means_vector_calculation(m_fr, weight_prop_basic)
baseline_fr = mean(temp_fr_abs[idx_ref])
temp_fr = temp_fr_abs .- baseline_fr

# -- CANADA --
m_ca = compute_general_climatology(data_folder_ca_31, weight_canada_31, YEARS, mode=:yearly)
temp_ca_abs = means_vector_calculation(m_ca, weight_canada_31)
baseline_ca = mean(temp_ca_abs[idx_ref])
temp_ca = temp_ca_abs .- baseline_ca


# 4. TENDANCES ET TESTS STATISTIQUES

# France
s_fr, i_fr, r2_fr, p_fr, se_fr, pred_fr = linear_trend(YEARS, temp_fr)
roll_fr = rolling_mean(temp_fr, 11)

# Canada
s_ca, i_ca, r2_ca, p_ca, se_ca, pred_ca = linear_trend(YEARS, temp_ca)
roll_ca = rolling_mean(temp_ca, 11)

# Conversion en °C/décennie
s_fr_dec = s_fr * 10
s_ca_dec = s_ca * 10

# Test statistique unilatéral (H1: La France chauffe plus vite que le Canada)
diff_slope = s_fr - s_ca
se_diff = sqrt(se_fr^2 + se_ca^2)
z_score = diff_slope / se_diff
p_diff = 1 - cdf(Normal(0, 1), z_score)


# 5. GRAPHIQUES ET SAUVEGARDE

println("Génération des graphiques...")

# France
p_fr_plot = plot(YEARS, temp_fr, label="Annuel", color=:blue, alpha=0.3, xlabel="Année", ylabel="Anomalie (°C)", title="France — Anomalies ERA5 (réf. 1960–1980)")
plot!(p_fr_plot, YEARS, roll_fr, label="Lissage 11 ans", color=:blue, lw=2)
plot!(p_fr_plot, YEARS, pred_fr, label=@sprintf("Tendance : %.2f °C/déc", s_fr_dec), color=:black, ls=:dash)
hline!(p_fr_plot, [0], color=:gray, ls=:dot, label="")
save_plot(p_fr_plot, joinpath(plot_dir, "france_tendance.png"))

# Canada
p_ca_plot = plot(YEARS, temp_ca, label="Annuel", color=:red, alpha=0.3, xlabel="Année", ylabel="Anomalie (°C)", title="Canada — Anomalies ERA5 (réf. 1960–1980)")
plot!(p_ca_plot, YEARS, roll_ca, label="Lissage 11 ans", color=:red, lw=2)
plot!(p_ca_plot, YEARS, pred_ca, label=@sprintf("Tendance : %.2f °C/déc", s_ca_dec), color=:black, ls=:dash)
hline!(p_ca_plot, [0], color=:gray, ls=:dot, label="")
save_plot(p_ca_plot, joinpath(plot_dir, "canada_tendance.png"))

# Comparaison
p_comp = plot(YEARS, roll_fr, label="France (lissé)", color=:blue, lw=2, xlabel="Année", ylabel="Anomalie (°C)", title="Comparaison France vs Canada — ERA5")
plot!(p_comp, YEARS, roll_ca, label="Canada (lissé)", color=:red, lw=2)
plot!(p_comp, YEARS, pred_fr, label=@sprintf("FR trend %.2f", s_fr_dec), color=:blue, ls=:dash)
plot!(p_comp, YEARS, pred_ca, label=@sprintf("CA trend %.2f", s_ca_dec), color=:red, ls=:dash)
hline!(p_comp, [0], color=:gray, ls=:dot, label="")
save_plot(p_comp, joinpath(plot_dir, "comparaison_tendances.png"))

# Panel Final 3x1
p_panel = plot(p_fr_plot, p_ca_plot, p_comp, layout=(3, 1), size=(1000, 1200))
save_plot(p_panel, joinpath(plot_dir, "panel_final_tendances.png"))

# 6. RÉSULTATS DANS LA CONSOLE

@printf("\n=== TENDANCES (réf. 1960–1980) ===\n")
@printf("France : %.3f °C/décennie (p=%.2e, R²=%.3f)\n", s_fr_dec, p_fr, r2_fr)
@printf("Canada : %.3f °C/décennie (p=%.2e, R²=%.3f)\n", s_ca_dec, p_ca, r2_ca)

ratio = s_fr / s_ca
@printf("\nLa France se réchauffe %.2f× plus vite que le Canada.\n", ratio)

@printf("\n=== TEST DIFFÉRENCE DE PENTES (FR > CA) ===\n")
@printf("Différence de pente (°C/an) : %.4f\n", diff_slope)
@printf("Erreur standard             : %.4f\n", se_diff)
@printf("z-score                     : %.2f\n", z_score)
@printf("p-valeur unilatérale        : p = %.3e\n", p_diff)

if p_diff < 0.05
    println("\n→ La différence est statistiquement significative : la France chauffe plus vite.")
else
    println("\n→ La différence n'est pas significative au seuil 5%.")
end