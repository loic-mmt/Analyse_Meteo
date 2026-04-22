# Vitesse de rechauffement ET COMPARAISON (FRANCE VS CANADA)

using Statistics
using Distributions
using Printf
using Plots
using StatsBase        # autocor()
using HypothesisTests  # JarqueBeraTest

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
# FONCTION 1 — Durbin-Watson
# Détecte l'autocorrélation des résidus
# DW ≈ 2 → OK   |   DW < 1.5 → autocorrélation positive (problème)
function durbin_watson(resid)
    return sum((resid[2:end] .- resid[1:end-1]).^2) / sum(resid.^2)
end

# FONCTION 2 — Correction de la variance pour l'autocorrélation
# Gonfle l'erreur standard si rho lag-1 est significatif
function effective_se(resid, se_ols, n)
    rho     = autocor(resid, [1])[1]
    n_eff   = max(n * (1 - rho) / (1 + rho), 2.0)
    inflation = sqrt(n / n_eff)
    return se_ols * inflation, rho, n_eff
end
# FONCTION 3 — Mann-Kendall (non-paramétrique)
# Robuste à l'autocorrélation, ne suppose pas la normalité
function mann_kendall(y)
    n = length(y)
    S = 0
    for i in 1:n-1
        for j in i+1:n
            S += sign(y[j] - y[i])
        end
    end
    varS  = n * (n - 1) * (2n + 5) / 18
    z     = S == 0 ? 0.0 : (S > 0 ? (S - 1) : (S + 1)) / sqrt(varS)
    p_bi  = 2 * (1 - cdf(Normal(0, 1), abs(z)))
    p_one = 1 - cdf(Normal(0, 1), z)
    return S, z, p_bi, p_one
end

# CALCUL DES RÉSIDUS 
resid_fr = temp_fr .- pred_fr
resid_ca = temp_ca .- pred_ca
n        = length(YEARS)
t95      = quantile(TDist(n - 2), 0.975)
 
# TEST 1 — Normalité des résidus (Jarque-Bera)
jb_fr   = JarqueBeraTest(resid_fr)
jb_ca   = JarqueBeraTest(resid_ca)
p_jb_fr = pvalue(jb_fr)
p_jb_ca = pvalue(jb_ca)
 
# TEST 2 — Autocorrélation (Durbin-Watson)
dw_fr = durbin_watson(resid_fr)
dw_ca = durbin_watson(resid_ca)

# CORRECTION — Erreurs standards corrigées de l'autocorrélation
se_fr_corr, rho_fr, neff_fr = effective_se(resid_fr, se_fr, n)
se_ca_corr, rho_ca, neff_ca = effective_se(resid_ca, se_ca, n)
 
ic_fr_corr = t95 * se_fr_corr * 10   # IC à 95% en °C/décennie
ic_ca_corr = t95 * se_ca_corr * 10
 
# TEST 3 — Mann-Kendall
S_fr, z_fr, p_mk_fr, p_mk_fr_one = mann_kendall(temp_fr)
S_ca, z_ca, p_mk_ca, p_mk_ca_one = mann_kendall(temp_ca)
 
# TEST 4 — Différence de pentes CORRIGÉE (SE effectives)
se_diff_corr = sqrt(se_fr_corr^2 + se_ca_corr^2)
z_diff_corr  = diff_slope / se_diff_corr
p_diff_corr  = 1 - cdf(Normal(0, 1), z_diff_corr)
 
# Chevauchement des IC
overlap = (s_fr_dec - ic_fr_corr < s_ca_dec + ic_ca_corr) &&
          (s_ca_dec - ic_ca_corr < s_fr_dec + ic_fr_corr)
 
# AFFICHAGE
@printf("TEST 1 — NORMALITÉ DES RÉSIDUS (Jarque-Bera) ===\n")
@printf("France : p = %.3f  → %s\n", p_jb_fr, p_jb_fr > 0.05 ? "✓ Normalité acceptée" : "✗ Normalité rejetée !")
@printf("Canada : p = %.3f  → %s\n", p_jb_ca, p_jb_ca > 0.05 ? "✓ Normalité acceptée" : "✗ Normalité rejetée !")
 
@printf(" TEST 2 — AUTOCORRÉLATION DES RÉSIDUS (Durbin-Watson) ===\n")
@printf("France : DW = %.3f  (rho lag-1 = %.3f, N_eff = %.1f / %d)\n", dw_fr, rho_fr, neff_fr, n)
@printf("Canada : DW = %.3f  (rho lag-1 = %.3f, N_eff = %.1f / %d)\n", dw_ca, rho_ca, neff_ca, n)
println("  DW ≈ 2 : pas d'autocorrélation  |  DW < 1.5 : autocorrélation détectée")
 
@printf(" INTERVALLES DE CONFIANCE 95%% CORRIGÉS ===\n")
@printf("France : %.3f ± %.3f °C/déc\n", s_fr_dec, ic_fr_corr)
@printf("Canada : %.3f ± %.3f °C/déc\n", s_ca_dec, ic_ca_corr)
println(overlap ?
    "   IC se chevauchent → différence non concluante sur ce critère" :
    "   IC ne se chevauchent pas")
 
@printf(" TEST 3 — MANN-KENDALL (non-paramétrique) ===\n")
@printf("France : S = %d, z = %.2f, p (unilatéral ↑) = %.4f  → %s\n",
    S_fr, z_fr, p_mk_fr_one, p_mk_fr_one < 0.05 ? "✓ Tendance significative" : "✗ Non significatif")
@printf("Canada : S = %d, z = %.2f, p (unilatéral ↑) = %.4f  → %s\n",
    S_ca, z_ca, p_mk_ca_one, p_mk_ca_one < 0.05 ? "✓ Tendance significative" : "✗ Non significatif")
 
@printf(" TEST 4 — DIFFÉRENCE DE PENTES CORRIGÉE ===\n")
@printf("OLS brut   : z = %.2f, p = %.4f\n", z_score, p_diff)
@printf("Corrigé AC : z = %.2f, p = %.4f\n", z_diff_corr, p_diff_corr)

crit1  = p_mk_fr_one < 0.05
crit2  = p_mk_ca_one < 0.05
crit3  = p_diff_corr < 0.05
crit4  = !overlap
n_crit = sum([crit1, crit2, crit3, crit4])
 
@printf("\n=== CONCLUSION RIGOUREUSE (%d/4 critères validés) ===\n", n_crit)
@printf("  [%s] Mann-Kendall France significatif\n",     crit1 ? "✓" : "✗")
@printf("  [%s] Mann-Kendall Canada significatif\n",     crit2 ? "✓" : "✗")
@printf("  [%s] Différence de pentes corrigée p<0.05\n", crit3 ? "✓" : "✗")
@printf("  [%s] IC 95%% ne se chevauchent pas\n",        crit4 ? "✓" : "✗")
 
if n_crit >= 3
    @printf("\n→ La France se réchauffe %.2f× plus vite que le Canada.\n", ratio)
    println("   Conclusion robuste (≥ 3 critères validés).")
elseif n_crit == 2
    println("\n→ Tendance France > Canada probable mais non fermement établie.")
else
    println("\n→ La différence France/Canada n'est PAS statistiquement établie.")
end
