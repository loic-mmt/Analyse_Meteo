using StatsPlots
using Statistics


# ================================================================================================
# =======
# Utilité fonctions
# =======
# ================================================================================================

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const plot_dir = joinpath(@__DIR__, "plot")
isdir(plot_dir) || mkpath(plot_dir)


# Chargement isolé des scripts pour éviter les collisions de noms de fonctions/constants.
module FranceCaniculeMod
include(joinpath(@__DIR__, "..", "France", "heatmap_canicule.jl"))
end

module FranceFroidMod
include(joinpath(@__DIR__, "..", "France", "heatmap_grand_froid_1950_2025.jl"))
end

module CanadaCaniculeMod
include(joinpath(@__DIR__, "..", "Canada", "heatmap_canicule.jl"))
end

module CanadaFroidMod
include(joinpath(@__DIR__, "..", "Canada", "heatmap_grand_froid_1950_2025_.jl"))
end


"""
    save_plot(plot_obj, plot_file)

Sauvegarde une figure en écrasant le fichier existant si besoin.
"""
function save_plot(plot_obj, plot_file::AbstractString)
    isfile(plot_file) && rm(plot_file)
    savefig(plot_obj, plot_file)
    return plot_file
end

"""
    moving_average_ignore_nan(values; window=7)

Calcule une moyenne glissante centrée en ignorant les NaN.
Renvoie un vecteur de même taille.
"""
function moving_average_ignore_nan(values::AbstractVector{<:Real}; window::Int=7)
    n = length(values)
    n == 0 && return Float64[]
    window = max(1, window)
    half_w = fld(window, 2)

    out = fill(NaN, n)
    arr = Float64.(values)

    for i in 1:n
        i1 = max(1, i - half_w)
        i2 = min(n, i + half_w)
        segment = @view arr[i1:i2]
        vals = segment[.!isnan.(segment)]
        out[i] = isempty(vals) ? NaN : mean(vals)
    end
    return out
end

"""
    align_year_series(years_a, vals_a, years_b, vals_b)

Aligne deux séries annuelles sur l'intersection des années.
"""
function align_year_series(
    years_a::AbstractVector{<:Integer},
    vals_a::AbstractVector{<:Real},
    years_b::AbstractVector{<:Integer},
    vals_b::AbstractVector{<:Real}
)
    map_a = Dict(Int(y) => Float64(v) for (y, v) in zip(years_a, vals_a))
    map_b = Dict(Int(y) => Float64(v) for (y, v) in zip(years_b, vals_b))

    common_years = sort(collect(intersect(Set(keys(map_a)), Set(keys(map_b)))))
    series_a = [map_a[y] for y in common_years]
    series_b = [map_b[y] for y in common_years]

    return common_years, series_a, series_b
end


# ================================================================================================
# Courbes comparées
# ================================================================================================

"""
    compute_canicule_comparison_series(; year_range=1950:2025, smooth_window=7)

Calcule les séries annuelles canicule Canada/France (jours moyens par pixel)
et leur moyenne lissée.
"""
function compute_canicule_comparison_series(; year_range=1950:2025, smooth_window::Int=7)
    _, fr_mean, fr_years = FranceCaniculeMod.compute_annual_canicule_days_maps_streaming(
        FranceCaniculeMod.data_folder_temporel,
        FranceCaniculeMod.weight_temporel,
        year_range;
        selected_months=collect(1:12),
        variable_name="t2m",
        day_threshold=36.0,
        night_threshold=21.0,
        min_consecutive_days=3
    )

    _, ca_mean, ca_years = CanadaCaniculeMod.compute_annual_canicule_days_maps_streaming(
        CanadaCaniculeMod.data_folder_temporel,
        CanadaCaniculeMod.weight_temporel,
        year_range;
        selected_months=collect(1:12),
        variable_name="t2m",
        min_consecutive_days=2,
        use_regional_thresholds=true
    )

    years, canada, france = align_year_series(ca_years, ca_mean, fr_years, fr_mean)
    canada_smooth = moving_average_ignore_nan(canada; window=smooth_window)
    france_smooth = moving_average_ignore_nan(france; window=smooth_window)

    return years, canada, france, canada_smooth, france_smooth
end

"""
    compute_cold_comparison_series(; year_range=1950:2025, smooth_window=7, cold_mode_france=:orange)

Calcule les séries annuelles grand froid Canada/France (jours moyens par pixel)
et leur moyenne lissée.
"""
function compute_cold_comparison_series(; year_range=1950:2025, smooth_window::Int=7, cold_mode_france::Symbol=:orange)
    _, fr_mean, fr_years = FranceFroidMod.compute_annual_cold_days_maps_streaming(
        FranceFroidMod.data_folder_temporel,
        FranceFroidMod.weight_temporel,
        year_range;
        cold_mode=cold_mode_france,
        selected_months=collect(1:12),
        variable_name="t2m",
        min_consecutive_days=1
    )

    _, ca_mean, ca_years = CanadaFroidMod.compute_annual_cold_days_maps_streaming(
        CanadaFroidMod.data_folder_temporel,
        CanadaFroidMod.weight_temporel,
        year_range;
        cold_mode=:regional,
        selected_months=collect(1:12),
        variable_name="t2m",
        min_consecutive_days=2
    )

    years, canada, france = align_year_series(ca_years, ca_mean, fr_years, fr_mean)
    canada_smooth = moving_average_ignore_nan(canada; window=smooth_window)
    france_smooth = moving_average_ignore_nan(france; window=smooth_window)

    return years, canada, france, canada_smooth, france_smooth
end

"""
    plot_canicule_comparison(...)

Trace la courbe comparée canicule :
- rouge : Canada
- bleu : France
- rouge pointillé : moyenne lissée Canada
- bleu pointillé : moyenne lissée France
"""
function plot_canicule_comparison(
    years::AbstractVector{<:Integer},
    canada::AbstractVector{<:Real},
    france::AbstractVector{<:Real},
    canada_smooth::AbstractVector{<:Real},
    france_smooth::AbstractVector{<:Real}
)
    p = plot(
        years,
        canada,
        color=:red,
        linewidth=2.2,
        label="Canada",
        xlabel="Année",
        ylabel="Jours de canicule / pixel",
        title="Comparaison Canada vs France - Canicule",
        legend=:topleft
    )
    plot!(p, years, france, color=:blue, linewidth=2.2, label="France")
    plot!(p, years, canada_smooth, color=:red, linestyle=:dash, linewidth=2.2, label="Moyenne lissée Canada")
    plot!(p, years, france_smooth, color=:blue, linestyle=:dash, linewidth=2.2, label="Moyenne lissée France")
    return p
end

"""
    plot_cold_comparison(...)

Trace la courbe comparée grand froid :
- rouge : Canada
- bleu : France
- rouge pointillé : moyenne lissée Canada
- bleu pointillé : moyenne lissée France
"""
function plot_cold_comparison(
    years::AbstractVector{<:Integer},
    canada::AbstractVector{<:Real},
    france::AbstractVector{<:Real},
    canada_smooth::AbstractVector{<:Real},
    france_smooth::AbstractVector{<:Real}
)
    p = plot(
        years,
        canada,
        color=:red,
        linewidth=2.2,
        label="Canada",
        xlabel="Année",
        ylabel="Jours de grand froid / pixel",
        title="Comparaison Canada vs France - Grand froid",
        legend=:topright
    )
    plot!(p, years, france, color=:blue, linewidth=2.2, label="France")
    plot!(p, years, canada_smooth, color=:red, linestyle=:dash, linewidth=2.2, label="Moyenne lissée Canada")
    plot!(p, years, france_smooth, color=:blue, linestyle=:dash, linewidth=2.2, label="Moyenne lissée France")
    return p
end

"""
    run_comparison_curves(; ...)

Pipeline complet :
1) calcule les séries canicule et grand froid (France + Canada),
2) génère les 2 graphiques comparatifs,
3) sauvegarde les figures.
"""
function run_comparison_curves(;
    year_range=1950:2025,
    smooth_window::Int=7,
    cold_mode_france::Symbol=:orange,
    canicule_output_file::String=joinpath(plot_dir, "comparaison_canicule_france_canada.png"),
    cold_output_file::String=joinpath(plot_dir, "comparaison_grand_froid_france_canada.png")
)
    println("Step 1: Canicule comparison data...")
    y_heat, ca_heat, fr_heat, ca_heat_s, fr_heat_s =
        compute_canicule_comparison_series(year_range=year_range, smooth_window=smooth_window)

    println("Step 2: Grand-froid comparison data...")
    y_cold, ca_cold, fr_cold, ca_cold_s, fr_cold_s =
        compute_cold_comparison_series(year_range=year_range, smooth_window=smooth_window, cold_mode_france=cold_mode_france)

    println("Step 3: Rendering plots...")
    p_heat = plot_canicule_comparison(y_heat, ca_heat, fr_heat, ca_heat_s, fr_heat_s)
    p_cold = plot_cold_comparison(y_cold, ca_cold, fr_cold, ca_cold_s, fr_cold_s)

    save_plot(p_heat, canicule_output_file)
    save_plot(p_cold, cold_output_file)

    println("Saved canicule comparison plot to: $canicule_output_file")
    println("Saved grand-froid comparison plot to: $cold_output_file")

    return p_heat, p_cold
end

"""
    main(args=ARGS)

Usage :
- `julia analyses_complementaires/Comparaison/courbes_comparees_canicule_grand_froid.jl`
- `julia analyses_complementaires/Comparaison/courbes_comparees_canicule_grand_froid.jl 1950 2025 7`
"""
function main(args=ARGS)
    start_year = length(args) >= 1 ? parse(Int, args[1]) : 1950
    end_year = length(args) >= 2 ? parse(Int, args[2]) : 2025
    smooth_window = length(args) >= 3 ? parse(Int, args[3]) : 7

    if end_year < start_year
        error("end_year ($end_year) doit être >= start_year ($start_year).")
    end

    run_comparison_curves(
        year_range=start_year:end_year,
        smooth_window=smooth_window
    )
    return nothing
end


# Équivalent Julia de `if __name__ == "__main__":`
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

