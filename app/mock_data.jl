module MockData

using Dates

import ..DataSource:
    AbstractTemperatureSource,
    TemperatureSeries,
    TemperatureField,
    available_years,
    default_months,
    temperature_field,
    temperature_series

export MockTemperatureSource, default_source

struct MockTemperatureSource <: AbstractTemperatureSource
    years::UnitRange{Int}
    lon::Vector{Float64}
    lat::Vector{Float64}
end

function MockTemperatureSource(years::UnitRange{Int}=1990:2024; grid::Tuple{Int, Int}=(60, 45))
    lon = collect(range(-5.0, stop=10.0, length=grid[1]))
    lat = collect(range(42.0, stop=51.0, length=grid[2]))
    return MockTemperatureSource(years, lon, lat)
end

default_source() = MockTemperatureSource()

available_years(source::MockTemperatureSource) = source.years

"""
    temperature_series(source; years, months)

Genere une serie de temperature France entiere pour les annees et mois demandes.
Les valeurs suivent une tendance de rechauffement douce avec un cycle annuel sinusoidal.
"""
function temperature_series(
    source::MockTemperatureSource;
    years::UnitRange{Int}=source.years,
    months::AbstractVector{<:Integer}=default_months(),
)::TemperatureSeries
    yspan = clamp_years(source, years)
    mspan = clamp_months(months)

    timestamps = Date[]
    values = Float64[]
    for y in yspan
        for m in mspan
            push!(timestamps, Date(y, m, 15))
            push!(values, synthetic_mean_temperature(source, y, m))
        end
    end

    return (timestamps=timestamps, values=values)
end

"""
    temperature_field(source; year, month)

Genere un champ 2D deterministe sur la grille France pour l'annee et le mois
demande. Un gradient nord-sud, une composante est-ouest et une tendance
saisonniere + long terme sont combines pour garder l'animation coherente.
"""
function temperature_field(
    source::MockTemperatureSource;
    year::Int,
    month::Int,
)::TemperatureField
    yr = clamp(year, first(source.years), last(source.years))
    mth = clamp(month, 1, 12)

    lon = source.lon
    lat = source.lat
    base = synthetic_mean_temperature(source, yr, mth)
    grid = [synthetic_cell_value(x, y, base, yr, mth) for x in lon, y in lat]

    return (year=yr, month=mth, grid=grid, extent=(lon=lon, lat=lat))
end

# -- aides ----------------------------------------------------------------------

function clamp_years(source::MockTemperatureSource, years::UnitRange{Int})
    start_year = clamp(first(years), first(source.years), last(source.years))
    stop_year = clamp(last(years), first(source.years), last(source.years))
    return start_year:stop_year
end

function clamp_months(months::AbstractVector{<:Integer})
    m = unique(round.(Int, months))
    return filter(x -> 1 <= x <= 12, sort(m))
end

function synthetic_mean_temperature(source::MockTemperatureSource, year::Int, month::Int)
    baseline = 11.5
    warming_trend = 0.035 * (year - first(source.years))
    seasonal = 7.5 * sin(2 * pi * (month - 1) / 12)
    intra_month = 0.4 * sin(0.2 * year + 0.8 * month)
    return baseline + warming_trend + seasonal + intra_month
end

function synthetic_cell_value(lon::Float64, lat::Float64, base::Float64, year::Int, month::Int)
    meridional_gradient = -0.45 * (lat - 46.5)
    zonal_gradient = 0.12 * (lon + 2.0)
    coastal_ridge = 0.6 * sin(0.6 * lon) * cos(0.4 * lat)
    seasonal_wave = 1.4 * sin(0.25 * lon + 0.15 * lat + month)
    slow_trend = 0.02 * (year - 1990)
    return base + meridional_gradient + zonal_gradient + coastal_ridge + seasonal_wave + slow_trend
end

end
