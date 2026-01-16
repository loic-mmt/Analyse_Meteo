module MeteoApp

using Dates
using Observables
using Bonito
using Bonito.DOM
using WGLMakie

include("data_source.jl")
include("mock_data.jl")

using .DataSource
using .MockData

export build_app, start_app, DEFAULT_PORT

const DEFAULT_PORT = 8050
const MONTH_LABELS_FR = ["janv", "fev", "mars", "avr", "mai", "juin", "juil", "aout", "sept", "oct", "nov", "dec"]

# -- aides ----------------------------------------------------------------------

decimal_year(d::Date) = year(d) + (dayofyear(d) - 1) / Dates.daysinyear(d)

function normalize_months(months::AbstractVector{<:Integer})
    cleaned = filter(x -> 1 <= x <= 12, unique(months))
    return isempty(cleaned) ? default_months() : sort(cleaned)
end

function range_from_slider(values::Vector{Int}, fallback::UnitRange{Int})
    if isempty(values)
        return fallback
    end
    lo = clamp(first(values), first(fallback), last(fallback))
    hi = clamp(length(values) >= 2 ? last(values) : lo, lo, last(fallback))
    return lo:hi
end

function month_labels_fr(months::AbstractVector{<:Integer})
    return [MONTH_LABELS_FR[m] for m in months if 1 <= m <= 12]
end

function averaged_field(source::AbstractTemperatureSource, year::Int, months::Vector{Int})
    months = normalize_months(months)
    frames = [temperature_field(source; year=year, month=m) for m in months]
    grids = map(f -> f.grid, frames)
    combined = reduce(+, grids) ./ length(grids)
    return (grid=combined, extent=first(frames).extent)
end

function build_figures(source::AbstractTemperatureSource)
    initial = temperature_field(source; year=first(available_years(source)), month=1)
    lon = initial.extent.lon
    lat = initial.extent.lat

    heat_obs = Observable(initial.grid)
    ts_x = Observable(Float64[])
    ts_y = Observable(Float64[])

    fig = Figure(resolution=(1180, 720))

    ax1 = Axis(fig[1, 1]; title="Temperature moyenne France (donnees fictives)", xlabel="Annee", ylabel="°C")
    lines!(ax1, ts_x, ts_y; color=:dodgerblue, linewidth=2)
    ax1.xtickformat = xs -> string.(round.(xs, digits=1))

    ax2 = Axis(fig[2, 1]; title="Carte de temperature France", xlabel="Longitude", ylabel="Latitude")
    hm = heatmap!(ax2, lon, lat, heat_obs; colormap=:thermal, interpolate=true)
    Colorbar(fig[2, 2], hm; label="°C")
    ax2.aspect = DataAspect()

    return fig, (ts_x=ts_x, ts_y=ts_y, heat=heat_obs, lon=lon, lat=lat, heatplot=hm)
end

function month_checkboxes()
    labels = MONTH_LABELS_FR
    boxes = [Bonito.Checkbox(true) for _ in 1:12]
    grid = DOM.div(
        (DOM.label(box, DOM.span(lbl)) for (box, lbl) in zip(boxes, labels))...;
        style="display:grid;grid-template-columns:repeat(4,1fr);gap:6px;font-size:12px;"
    )
    return boxes, grid
end

function selected_months(boxes)
    return [i for (i, box) in enumerate(boxes) if box.value[]]
end

function ensure_checked!(boxes, months)
    active = Set(months)
    for (i, box) in enumerate(boxes)
        box.value[] = i in active
    end
end

# -- app composition ------------------------------------------------------------

function build_app(; source::AbstractTemperatureSource=default_source())
    WGLMakie.activate!()

    return App() do session
        fig, handles = build_figures(source)

        year_span = available_years(source)
        year_slider = Bonito.RangeSlider(year_span; attributes=Dict(:step => 1))
        year_slider.range[] = year_span
        year_slider.value[] = [first(year_span), last(year_span)]

        frame_slider = Bonito.Slider(collect(year_span); attributes=Dict(:step => 1))
        frame_slider.value[] = first(year_span)
        play_button = Bonito.Button("Lecture / Pause")
        playing = Observable(false)
        frame_year = Observable(first(year_span))

        boxes, boxes_dom = month_checkboxes()
        selected_months_obs = Observable(default_months())

        series_info = Observable("")
        map_info = Observable("")

        function refresh_series()
            yrange = range_from_slider(year_slider.value[], year_span)
            months = normalize_months(selected_months_obs[])
            data = temperature_series(source; years=yrange, months=months)
            handles.ts_x[] = [decimal_year(d) for d in data.timestamps]
            handles.ts_y[] = data.values
            month_label = join(month_labels_fr(months), ", ")
            series_info[] = "Annees $(first(yrange))-$(last(yrange)), mois : $(month_label)"
        end

        function refresh_map()
            months = normalize_months(selected_months_obs[])
            frame = averaged_field(source, frame_year[], months)
            handles.heat[] = frame.grid
            month_label = join(month_labels_fr(months), ", ")
            map_info[] = "Carte annee $(frame_year[]), mois : $(month_label)"
        end

        refresh_series()
        refresh_map()

        on(year_slider.value) do vals
            yrange = range_from_slider(vals, year_span)
            frame_slider.values[] = collect(yrange)
            frame_year[] = clamp(frame_year[], first(yrange), last(yrange))
            refresh_series()
            refresh_map()
        end

        on(frame_slider.value) do val
            frame_year[] = val
        end

        on(frame_year) do _
            refresh_map()
        end

        for cb in boxes
            on(cb.value) do _
                months = selected_months(boxes)
                if isempty(months)
                    months = default_months()
                    ensure_checked!(boxes, months)
                end
                selected_months_obs[] = months
                refresh_series()
                refresh_map()
            end
        end

        on(play_button) do _
            playing[] = !playing[]
        end

        @async begin
            while true
                if playing[]
                    yrange = range_from_slider(year_slider.value[], year_span)
                    next_year = frame_year[] >= last(yrange) ? first(yrange) : frame_year[] + 1
                    frame_slider.value[] = next_year
                    frame_year[] = next_year
                end
                sleep(0.9)
            end
        end

        control_panel = DOM.div(
            DOM.h2("Tableau de bord temperature France (donnees fictives)"),
            DOM.div(DOM.strong("Plage d'annees"), year_slider; class="control-row"),
            DOM.div(DOM.strong("Mois"), boxes_dom; class="control-row"),
            DOM.div(
                DOM.strong("Annee pour la carte"),
                DOM.span("Utilisez le curseur ou Lecture pour animer"),
                DOM.div(frame_slider, play_button, DOM.span(frame_year); class="frame-row");
                class="control-row"
            ),
            DOM.div(DOM.span(series_info); class="info-row"),
            DOM.div(DOM.span(map_info); class="info-row");
            class="controls"
        )

        layout = DOM.div(
            DOM.style("""
            body { font-family: Arial, sans-serif; background: #f8f9fb; }
            .controls { display:flex; flex-direction:column; gap:12px; padding:12px; background:white; border:1px solid #e0e4ea; border-radius:8px; }
            .control-row { display:flex; flex-direction:column; gap:4px; }
            .frame-row { display:flex; align-items:center; gap:10px; flex-wrap:wrap; }
            .info-row { color:#4a5568; font-size:13px; }
            .main { display:grid; grid-template-columns: 360px 1fr; gap:16px; align-items:start; }
            .plot-card { background:white; border:1px solid #e0e4ea; border-radius:8px; padding:8px; }
            @media (max-width: 900px) {
                .main { grid-template-columns: 1fr; }
            }
            """),
            DOM.div(
                control_panel,
                DOM.div(fig; class="plot-card"),
                class="main"
            )
        )

        return layout
    end
end

function start_app(; source::AbstractTemperatureSource=default_source(), host::AbstractString="127.0.0.1", port::Integer=DEFAULT_PORT)
    app = build_app(; source=source)
    server = Server(app, host, port)
    Bonito.HTTPServer.start(server)
    println("Application en cours d'execution sur http://$(server.url):$(server.port)/")
    return (app=app, server=server)
end

if abspath(PROGRAM_FILE) == @__FILE__
    start_app()
    wait()
end

end
