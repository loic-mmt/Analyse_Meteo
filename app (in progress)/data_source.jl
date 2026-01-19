module DataSource

using Dates

export AbstractTemperatureSource,
    TemperatureSeries,
    TemperatureField,
    temperature_series,
    temperature_field,
    available_years,
    default_months

# Definition minimale d'interface pour toute source de temperatures que l'appli peut appeler.
abstract type AbstractTemperatureSource end

const TemperatureSeries = NamedTuple{(:timestamps, :values), Tuple{Vector{Date}, Vector{Float64}}}
const TemperatureField = NamedTuple{(:year, :month, :grid, :extent), Tuple{Int, Int, Matrix{Float64}, NamedTuple{(:lon, :lat), Tuple{Vector{Float64}, Vector{Float64}}}}}

"""
    default_months()

Renvoie la liste par defaut des mois a tracer (janvier a decembre).
"""
default_months() = collect(1:12)

"""
    available_years(source)

Renvoie l'intervalle d'annees valides expose par la source.
"""
function available_years(::AbstractTemperatureSource)::UnitRange{Int}
    error("available_years is not implemented for this data source")
end

"""
    temperature_series(source; years, months) -> TemperatureSeries

Recupere une serie temporelle couvrant les annees et mois demandes. Le resultat
est un `NamedTuple` avec les horodatages (Vector{Date}) et les valeurs (Vector{Float64}).
"""
function temperature_series(
    ::AbstractTemperatureSource;
    years::UnitRange{Int},
    months::AbstractVector{<:Integer},
)::TemperatureSeries
    error("temperature_series is not implemented for this data source")
end

"""
    temperature_field(source; year, month) -> TemperatureField

Recupere un raster de temperature pour l'`year` et le `month` demandes. Le retour
est un `NamedTuple` avec l'annee, le mois, la grille (matrice) et les vecteurs
d'axes lon/lat.
"""
function temperature_field(
    ::AbstractTemperatureSource;
    year::Int,
    month::Int,
)::TemperatureField
    error("temperature_field is not implemented for this data source")
end

end
