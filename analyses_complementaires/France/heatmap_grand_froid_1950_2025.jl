using NCDatasets # Chargement des fichiers .nc
using Glob # Recherche de fichiers par motif
using StatsPlots # heatmap / savefig
using Dates # day/hour sur DateTime
using ArchGDAL # Lecture shapefile
using Proj # Reprojection vers WGS84


# ================================================================================================
# =======
# Utilité fonctions
# =======
# ================================================================================================

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# Données température
const data_folder_precise = joinpath(PROJECT_ROOT, "data", "france", "france_9")
const data_folder_basic = joinpath(PROJECT_ROOT, "data", "france", "france_9_31")
const data_folder_temporel = joinpath(PROJECT_ROOT, "data", "raw_monthly_combined", "precise")
# Fichiers de poids
const weight_prop_basic = joinpath(PROJECT_ROOT, "data", "masks", "weights_france_31.nc")
const weight_prop_precise = joinpath(PROJECT_ROOT, "data", "masks", "weights_france_9.nc")
const weight_temporel = joinpath(PROJECT_ROOT, "data", "masks", "weights_prop_precise.nc")

# Dossier de sortie des figures
const plot_dir = joinpath(@__DIR__, "plot")
isdir(plot_dir) || mkpath(plot_dir)

# Shapefile du contour France
const FRANCE_SHP = joinpath(PROJECT_ROOT, "data", "shapefiles", "region.shp")
const OUTLINE_CACHE = Dict{String, Vector{Tuple{Vector{Float64}, Vector{Float64}}}}()


"""
    save_plot(plot_obj, plot_file)

Sauvegarde une figure en écrasant le fichier existant si besoin.
"""
function save_plot(plot_obj, plot_file)
    isfile(plot_file) && rm(plot_file)
    savefig(plot_obj, plot_file)
end

"""
    print_progress(prefix, current, total; width=30)

Affiche une barre de progression simple en ligne (`\\r`) pour suivre
l'avancement des boucles longues.
"""
function print_progress(prefix::AbstractString, current::Int, total::Int; width::Int=30)
    total <= 0 && return
    ratio = clamp(current / total, 0.0, 1.0)
    filled = round(Int, ratio * width)
    bar = repeat("#", filled) * repeat("-", width - filled)
    pct = round(ratio * 100; digits=1)
    print("\r$(prefix) [$(bar)] $(pct)% ($(current)/$(total))")
    flush(stdout)
    if current == total
        println()
    end
end

"""
    compute_day_stats!(day_cube, day_max, day_min, valid_day; kelvin_to_celsius=true)

Calcule en place Tmax/Tmin journaliers et un masque de validité pixel, avec
conversion optionnelle Kelvin -> Celsius. Conçu pour limiter les allocations.
"""
function compute_day_stats!(
    day_cube,
    day_max::Matrix{Float64},
    day_min::Matrix{Float64},
    valid_day::BitMatrix;
    kelvin_to_celsius::Bool=true
)
    fill!(day_max, -Inf)
    fill!(day_min, Inf)
    fill!(valid_day, true)

    @inbounds for h in axes(day_cube, 3)
        frame = @view day_cube[:, :, h]
        for j in axes(frame, 2)
            for i in axes(frame, 1)
                v0 = frame[i, j]
                if ismissing(v0)
                    valid_day[i, j] = false
                    continue
                end
                v = Float64(v0)
                kelvin_to_celsius && (v -= 273.15)
                if isnan(v)
                    valid_day[i, j] = false
                    continue
                end
                if v > day_max[i, j]
                    day_max[i, j] = v
                end
                if v < day_min[i, j]
                    day_min[i, j] = v
                end
            end
        end
    end
    return nothing
end


# ================================================================================================
# Contours et masque France
# ================================================================================================

"""
    country_mask_for_map(weights, target_size)

Construit un masque binaire `Bool` (zone valide / zone hors carte) aligné
avec la taille de la matrice à tracer.
"""
function country_mask_for_map(weights::AbstractMatrix{<:Real}, target_size::Tuple{Int, Int})
    mask = weights .> 0.0
    if size(mask) == target_size
        return mask
    end

    mask_t = permutedims(mask)
    if size(mask_t) == target_size
        return mask_t
    end

    error("Incompatible sizes between weights $(size(weights)) and map $(target_size).")
end


"""
    shapefile_to_wgs84(srs)

Crée la transformation de projection vers EPSG:4326 (lon/lat),
format attendu pour superposer le contour sur la heatmap.
"""
function shapefile_to_wgs84(srs)
    crs_wkt = try
        ArchGDAL.toWKT(srs)
    catch
        ""
    end
    crs_src = !isempty(strip(crs_wkt)) ? crs_wkt : (try
        ArchGDAL.toPROJ4(srs)
    catch
        ""
    end)

    isempty(strip(crs_src)) && error("Unable to export shapefile CRS to WKT/PROJ4.")
    return Proj.Transformation(crs_src, "EPSG:4326", always_xy=true)
end


"""
    collect_outline_segments!(segs, geom, pj)

Parcourt récursivement une géométrie et ajoute les segments (lon/lat)
au vecteur `segs`.
"""
function collect_outline_segments!(segs::Vector{Tuple{Vector{Float64}, Vector{Float64}}}, geom, pj)
    gname = uppercase(ArchGDAL.geomname(geom))

    if startswith(gname, "LINESTRING") && !occursin("MULTI", gname)
        n = ArchGDAL.ngeom(geom)
        if n >= 2
            max_points = 2000
            step = max(1, cld(n, max_points))
            idxs = collect(0:step:(n - 1))
            (last(idxs) == n - 1) || push!(idxs, n - 1)

            xs = Vector{Float64}(undef, length(idxs))
            ys = Vector{Float64}(undef, length(idxs))
            for (k, i) in enumerate(idxs)
                x, y, _ = ArchGDAL.getpoint(geom, i)
                lon, lat = pj(Float64(x), Float64(y))
                xs[k] = Float64(lon)
                ys[k] = Float64(lat)
            end
            push!(segs, (xs, ys))
        end
        return
    end

    for i in 0:(ArchGDAL.ngeom(geom) - 1)
        collect_outline_segments!(segs, ArchGDAL.getgeom(geom, i), pj)
    end
end


@inline function segment_length(seg::Tuple{Vector{Float64}, Vector{Float64}})
    xs, ys = seg
    acc = 0.0
    for i in 2:length(xs)
        acc += hypot(xs[i] - xs[i - 1], ys[i] - ys[i - 1])
    end
    return acc
end


"""
    load_outline_segments(shp_path)

Charge les segments du contour France depuis un shapefile.
Le résultat est mis en cache pour accélérer les appels suivants.
"""
function load_outline_segments(shp_path::AbstractString)
    key = String(shp_path)
    if haskey(OUTLINE_CACHE, key)
        return OUTLINE_CACHE[key]
    end

    isfile(shp_path) || error("Shapefile not found: $shp_path")

    segments = ArchGDAL.read(shp_path) do dataset
        layer = ArchGDAL.getlayer(dataset, 0)
        srs = ArchGDAL.getspatialref(layer)
        pj = shapefile_to_wgs84(srs)

        geoms = ArchGDAL.IGeometry[]
        for feat in layer
            push!(geoms, ArchGDAL.clone(ArchGDAL.getgeom(feat)))
        end
        isempty(geoms) && error("No geometries found in $shp_path")

        u = geoms[1]
        for k in 2:length(geoms)
            u = ArchGDAL.union(u, geoms[k])
        end

        boundary = ArchGDAL.boundary(u)
        segs = Tuple{Vector{Float64}, Vector{Float64}}[]
        collect_outline_segments!(segs, boundary, pj)

        # On garde les segments significatifs pour limiter le coût de tracé.
        lengths = [segment_length(seg) for seg in segs]
        keep = findall(>(0.05), lengths)
        if !isempty(keep)
            segs = segs[keep]
            lengths = lengths[keep]
        end

        max_segments = 1200
        if length(segs) > max_segments
            order = partialsortperm(lengths, 1:max_segments; rev=true)
            segs = segs[order]
        end

        segs
    end

    OUTLINE_CACHE[key] = segments
    return segments
end


"""
    add_country_outline!(p, outline_segments; linecolor=:black, linewidth=1.2)

Ajoute le contour du pays sur une figure existante.
"""
function add_country_outline!(p, outline_segments::Vector{Tuple{Vector{Float64}, Vector{Float64}}}; linecolor=:black, linewidth=1.2)
    for (xs, ys) in outline_segments
        plot!(p, xs, ys; color=linecolor, linewidth=linewidth, label=false)
    end
    return p
end


"""
    prepare_heatmap_axes(lons, lats, z)

Trie les axes lon/lat et retourne les index à appliquer à `z` avant affichage.
"""
function prepare_heatmap_axes(lons::AbstractVector, lats::AbstractVector, z::AbstractMatrix)
    nrows, ncols = size(z)
    nrows == length(lats) || error("Heatmap rows ($nrows) != latitude length ($(length(lats))).")
    ncols == length(lons) || error("Heatmap cols ($ncols) != longitude length ($(length(lons))).")

    lon_idx = sortperm(lons)
    lat_idx = sortperm(lats)
    x_plot = Float64.(lons[lon_idx])
    y_plot = Float64.(lats[lat_idx])
    return x_plot, y_plot, lon_idx, lat_idx
end


# ================================================================================================
# Préparer les données
# ================================================================================================

"""
    compute_general_climatology(data_folder, weights_file, year_range; ...)

Orchestrateur de préparation des données ERA5 :
1. accumulation depuis les NetCDF horaires,
2. calcul des moyennes,
3. conversion Kelvin -> Celsius,
4. application du masque spatial.

Le mode attendu pour le grand froid est `mode=:hourly`.

# Retourne
- `final_cube::Array{Float64,3}` : cube `[lon, lat, time]` en Celsius.
"""
function compute_general_climatology(
    data_folder::String,
    weights_file::String,
    year_range;
    mode::Symbol=:hourly,
    selected_months=collect(1:12),
    selected_days=nothing,
    selected_hours=0:23,
    variable_name="t2m"
)
    if year_range isa Integer
        year_range = [year_range]
    end
    if selected_months isa Integer
        selected_months = [selected_months]
    end
    if selected_hours isa Integer
        selected_hours = [selected_hours]
    end
    if selected_days isa Integer
        selected_days = [selected_days]
    end

    ds = NCDataset(weights_file)
    weights = ds["final_weights"][:, :]
    close(ds)

    println("Step 1: Accumulating data...")
    sums, counts, _lons, _lats = accumulate_data(
        data_folder, year_range, mode, selected_months, selected_days, selected_hours, variable_name
    )

    println("Step 2: Computing means...")
    final_cube, _valid_times = finalize_cube(sums, counts, weights, mode)

    return final_cube
end


"""
    accumulate_data(data_folder, year_range, mode, months, days, hours, var_name)

Lit les fichiers NetCDF sans tout charger en mémoire.
Accumulates `sum` et `count` par clé temporelle selon le `mode`.
"""
function accumulate_data(data_folder, year_range, mode, months, days, hours, var_name)
    sums_dict = Dict{Any, Matrix{Float64}}()
    counts_dict = Dict{Any, Matrix{Int}}()
    lons, lats = nothing, nothing
    total_steps = length(year_range) * length(months)
    step = 0

    for year in year_range
        for month in months
            step += 1
            month_str = lpad(month, 2, '0')
            files = sort(glob("*$(year)_$(month_str)*.nc", data_folder))
            if isempty(files)
                print_progress("Accumulation", step, total_steps)
                continue
            end

            for file in files
                NCDataset(file) do ds
                    if isnothing(lons)
                        lons, lats = ds["longitude"][:], ds["latitude"][:]
                    end

                    times = ds["valid_time"][:]
                    indices = findall(t ->
                        (isnothing(days) || day(t) in days) &&
                        (hour(t) in hours),
                        times
                    )

                    if !isempty(indices)
                        data_slice = ds[var_name][:, :, indices]

                        for t in 1:length(indices)
                            time_val = times[indices[t]]
                            key = get_binning_key(mode, year, month, time_val)

                            if !haskey(sums_dict, key)
                                dims = size(data_slice)[1:2]
                                sums_dict[key] = zeros(Float64, dims)
                                counts_dict[key] = zeros(Int, dims)
                            end

                            frame = data_slice[:, :, t]
                            valid = map(x -> !ismissing(x) && !isnan(Float64(x)), frame)

                            if any(valid)
                                sums_dict[key][valid] .+= Float64.(frame[valid])
                                counts_dict[key][valid] .+= 1
                            end
                        end
                    end
                end
            end
            print_progress("Accumulation", step, total_steps)
        end
    end

    return sums_dict, counts_dict, lons, lats
end


"""
    get_binning_key(mode, year, month, time_val)

Construit la clé temporelle d'agrégation.
"""
function get_binning_key(mode, year, month, time_val)
    if mode == :monthly
        return Date(year, month, 1)
    elseif mode == :daily
        return Date(year, month, day(time_val))
    elseif mode == :hourly
        return time_val
    elseif mode == :yearly
        return year
    else
        return "Total"
    end
end


"""
    finalize_cube(sums_dict, counts_dict, weights, mode)

Transforme les accumulations en cube 3D en Celsius avec masque spatial.
"""
function finalize_cube(sums_dict, counts_dict, weights, mode)
    all_keys = collect(keys(sums_dict))
    mode != :total && sort!(all_keys)

    valid_keys = filter(k -> any(counts_dict[k] .> 0), all_keys)
    if isempty(valid_keys)
        println("Warning: No valid data found.")
        return Array{Float64}(undef, 0, 0, 0), []
    end

    first_data = sums_dict[valid_keys[1]]
    n_lon, n_lat = size(first_data)

    weight_map =
        if size(weights) == (n_lon, n_lat)
            weights
        elseif size(weights) == (n_lat, n_lon)
            permutedims(weights)
        else
            error("Incompatible sizes between weights $(size(weights)) and data ($((n_lon, n_lat))).")
        end

    visual_mask = fill(NaN, size(weight_map))
    visual_mask[weight_map .> 0.0] .= 1.0

    n_time = length(valid_keys)
    final_3d = fill(NaN, n_lon, n_lat, n_time)

    for (i, k) in enumerate(valid_keys)
        sum_grid = sums_dict[k]
        count_grid = counts_dict[k]

        mean_grid = fill(NaN, n_lon, n_lat)
        valid = count_grid .> 0
        mean_grid[valid] = sum_grid[valid] ./ count_grid[valid]

        # Kelvin -> Celsius
        mean_grid .-= 273.15

        # Masque spatial France
        mean_grid .*= visual_mask
        final_3d[:, :, i] = mean_grid
    end

    return final_3d, valid_keys
end


# ================================================================================================
# Détection grand froid
# ================================================================================================

"""
    normalize_cold_mode(mode)

Normalise le mode d'alerte grand froid vers `:jaune`, `:orange` ou `:rouge`.
"""
function normalize_cold_mode(mode::Symbol)
    m = Symbol(lowercase(String(mode)))
    if m in (:jaune, :yellow, :niveau1, :level1, :n1)
        return :jaune
    elseif m in (:orange, :niveau2, :level2, :n2)
        return :orange
    elseif m in (:rouge, :red, :niveau3, :level3, :n3)
        return :rouge
    else
        error("Mode invalide: $mode. Valeurs attendues: :jaune, :orange, :rouge")
    end
end

"""
    cold_mode_label(mode)

Renvoie le libellé texte du mode d'alerte.
"""
cold_mode_label(mode::Symbol) = String(normalize_cold_mode(mode))

"""
    is_cold_day(tmax, tmin, cold_mode)

Qualification journalière du grand froid selon le mode choisi.
"""
function is_cold_day(tmax::Float64, tmin::Float64, cold_mode::Symbol)
    mode = normalize_cold_mode(cold_mode)
    if mode == :jaune
        # Entre 0°C et -5°C, avec possibilité de descendre jusqu'à -10°C la nuit.
        return (tmax <= 0.0) && (tmin <= -5.0) && (tmin > -10.0)
    elseif mode == :orange
        # Entre -10°C et -18°C.
        return (tmin <= -10.0) && (tmin > -18.0)
    else
        # Sous -18°C.
        return tmin <= -18.0
    end
end

"""
    calculate_cold_days(data_3d, weights_file; ...)

Calcule le total d'heures de grand froid par pixel selon le mode choisi.
"""
function calculate_cold_days(
    data_3d::AbstractArray{<:Union{Missing, Real}, 3},
    weights_file;
    cold_mode::Symbol=:orange,
    min_consecutive_days::Int=1,
    mask_threshold::Float64=0.0,
    block_hours::Int=24
)
    mode = normalize_cold_mode(cold_mode)

    ds = NCDataset(weights_file)
    weights = ds["weights_frac"][:, :]
    close(ds)

    n_lon, n_lat, n_time = size(data_3d)
    n_time % block_hours == 0 || error("La dimension temporelle ($n_time) n'est pas divisible par $block_hours.")
    n_days = div(n_time, block_hours)

    weight_mask =
        if size(weights) == (n_lat, n_lon)
            weights
        elseif size(weights) == (n_lon, n_lat)
            permutedims(weights)
        else
            error("Dimensions incompatibles entre poids $(size(weights)) et data_3d ($(n_lon), $(n_lat), $n_time).")
        end

    cold_hours_grid = fill(NaN, n_lon, n_lat)

    println("Starting grand-froid calculations (mode: $(cold_mode_label(mode)))...")

    for i in 1:n_lon
        for j in 1:n_lat
            if weight_mask[j, i] <= mask_threshold
                continue
            end

            y_temps = @view data_3d[i, j, :]
            is_cold = falses(n_days)

            for d in 1:n_days
                i1 = block_hours * (d - 1) + 1
                i2 = i1 + (block_hours - 1)
                vals_day = @view y_temps[i1:i2]

                if any(ismissing, vals_day) || any(x -> isnan(Float64(x)), vals_day)
                    continue
                end

                vals_float = Float64.(vals_day)
                tmax = maximum(vals_float)
                tmin = minimum(vals_float)
                is_cold[d] = is_cold_day(tmax, tmin, mode)
            end

            run_len = 0
            cold_days = 0
            for d in eachindex(is_cold)
                if is_cold[d]
                    run_len += 1
                else
                    if run_len >= min_consecutive_days
                        cold_days += run_len
                    end
                    run_len = 0
                end
            end
            if run_len >= min_consecutive_days
                cold_days += run_len
            end

            cold_hours_grid[i, j] = cold_days * block_hours
        end
        print_progress("Grand froid", i, n_lon)
    end

    return cold_hours_grid
end

"""
    finalize_open_runs!(run_len, event_days, min_consecutive_days)

Ajoute les séquences encore ouvertes à la fin d'un bloc temporel.
"""
function finalize_open_runs!(run_len::AbstractMatrix{<:Integer}, event_days::AbstractMatrix{<:Integer}, min_consecutive_days::Int)
    @inbounds for j in axes(run_len, 2)
        for i in axes(run_len, 1)
            rl = run_len[i, j]
            if rl >= min_consecutive_days
                event_days[i, j] += rl
            end
        end
    end
    return nothing
end

"""
    calculate_cold_days_streaming(data_folder, weights_file, year_range; ...)

Version basse mémoire : lit les NetCDF mois par mois et calcule les heures
de grand froid sans construire de cube 3D global.
"""
function calculate_cold_days_streaming(
    data_folder::String,
    weights_file::String,
    year_range;
    cold_mode::Symbol=:orange,
    selected_months=collect(1:12),
    variable_name::String="t2m",
    min_consecutive_days::Int=1,
    mask_threshold::Float64=0.0,
    day_hours::Int=24
)
    mode = normalize_cold_mode(cold_mode)

    if year_range isa Integer
        year_range = [year_range]
    end
    if selected_months isa Integer
        selected_months = [selected_months]
    end

    first_file = nothing
    for year in year_range, month in selected_months
        month_str = lpad(month, 2, '0')
        files = sort(glob("*$(year)_$(month_str)*.nc", data_folder))
        if !isempty(files)
            first_file = files[1]
            break
        end
    end
    isnothing(first_file) && error("Aucun fichier NetCDF trouvé dans $data_folder pour la période demandée.")

    n_lon, n_lat = NCDataset(first_file) do ds0
        dims = size(ds0[variable_name])
        dims[1], dims[2]
    end

    dsw = NCDataset(weights_file)
    weights = dsw["weights_frac"][:, :]
    close(dsw)
    land_mask =
        if size(weights) == (n_lon, n_lat)
            weights .> mask_threshold
        elseif size(weights) == (n_lat, n_lon)
            permutedims(weights) .> mask_threshold
        else
            error("Dimensions incompatibles entre poids $(size(weights)) et grille ($((n_lon, n_lat))).")
        end

    run_len = zeros(Int16, n_lon, n_lat)
    cold_days_count = zeros(Int32, n_lon, n_lat)
    day_max = zeros(Float64, n_lon, n_lat)
    day_min = zeros(Float64, n_lon, n_lat)
    valid_day = falses(n_lon, n_lat)

    total_steps = length(year_range) * length(selected_months)
    step = 0
    println("Starting streaming grand-froid calculations (mode: $(cold_mode_label(mode)))...")

    for year in year_range
        for month in selected_months
            step += 1
            month_str = lpad(month, 2, '0')
            files = sort(glob("*$(year)_$(month_str)*.nc", data_folder))

            if isempty(files)
                finalize_open_runs!(run_len, cold_days_count, min_consecutive_days)
                run_len .= 0
                print_progress("Grand froid stream", step, total_steps)
                continue
            end

            for file in files
                NCDataset(file) do ds
                    n_t = size(ds[variable_name], 3)
                    n_full_days = fld(n_t, day_hours)
                    day_start = 1

                    if mode == :jaune
                        for _ in 1:n_full_days
                            day_end = day_start + day_hours - 1
                            idx = day_start:day_end
                            day_cube = ds[variable_name][:, :, idx]
                            compute_day_stats!(day_cube, day_max, day_min, valid_day; kelvin_to_celsius=true)

                            @inbounds for j in 1:n_lat
                                for i in 1:n_lon
                                    event = land_mask[i, j] && valid_day[i, j] &&
                                            (day_max[i, j] <= 0.0) &&
                                            (day_min[i, j] <= -5.0) &&
                                            (day_min[i, j] > -10.0)
                                    if event
                                        run_len[i, j] += 1
                                    else
                                        if run_len[i, j] >= min_consecutive_days
                                            cold_days_count[i, j] += run_len[i, j]
                                        end
                                        run_len[i, j] = 0
                                    end
                                end
                            end

                            day_start = day_end + 1
                        end
                    elseif mode == :orange
                        for _ in 1:n_full_days
                            day_end = day_start + day_hours - 1
                            idx = day_start:day_end
                            day_cube = ds[variable_name][:, :, idx]
                            compute_day_stats!(day_cube, day_max, day_min, valid_day; kelvin_to_celsius=true)

                            @inbounds for j in 1:n_lat
                                for i in 1:n_lon
                                    event = land_mask[i, j] && valid_day[i, j] &&
                                            (day_min[i, j] <= -10.0) &&
                                            (day_min[i, j] > -18.0)
                                    if event
                                        run_len[i, j] += 1
                                    else
                                        if run_len[i, j] >= min_consecutive_days
                                            cold_days_count[i, j] += run_len[i, j]
                                        end
                                        run_len[i, j] = 0
                                    end
                                end
                            end

                            day_start = day_end + 1
                        end
                    else
                        for _ in 1:n_full_days
                            day_end = day_start + day_hours - 1
                            idx = day_start:day_end
                            day_cube = ds[variable_name][:, :, idx]
                            compute_day_stats!(day_cube, day_max, day_min, valid_day; kelvin_to_celsius=true)

                            @inbounds for j in 1:n_lat
                                for i in 1:n_lon
                                    event = land_mask[i, j] && valid_day[i, j] &&
                                            (day_min[i, j] <= -18.0)
                                    if event
                                        run_len[i, j] += 1
                                    else
                                        if run_len[i, j] >= min_consecutive_days
                                            cold_days_count[i, j] += run_len[i, j]
                                        end
                                        run_len[i, j] = 0
                                    end
                                end
                            end

                            day_start = day_end + 1
                        end
                    end

                    # Fallback robuste si un fichier contient une journée incomplète en fin.
                    if day_start <= n_t
                        idx = day_start:n_t
                        day_cube = ds[variable_name][:, :, idx]
                        compute_day_stats!(day_cube, day_max, day_min, valid_day; kelvin_to_celsius=true)

                        if mode == :jaune
                            @inbounds for j in 1:n_lat
                                for i in 1:n_lon
                                    event = land_mask[i, j] && valid_day[i, j] &&
                                            (day_max[i, j] <= 0.0) &&
                                            (day_min[i, j] <= -5.0) &&
                                            (day_min[i, j] > -10.0)
                                    if event
                                        run_len[i, j] += 1
                                    else
                                        if run_len[i, j] >= min_consecutive_days
                                            cold_days_count[i, j] += run_len[i, j]
                                        end
                                        run_len[i, j] = 0
                                    end
                                end
                            end
                        elseif mode == :orange
                            @inbounds for j in 1:n_lat
                                for i in 1:n_lon
                                    event = land_mask[i, j] && valid_day[i, j] &&
                                            (day_min[i, j] <= -10.0) &&
                                            (day_min[i, j] > -18.0)
                                    if event
                                        run_len[i, j] += 1
                                    else
                                        if run_len[i, j] >= min_consecutive_days
                                            cold_days_count[i, j] += run_len[i, j]
                                        end
                                        run_len[i, j] = 0
                                    end
                                end
                            end
                        else
                            @inbounds for j in 1:n_lat
                                for i in 1:n_lon
                                    event = land_mask[i, j] && valid_day[i, j] &&
                                            (day_min[i, j] <= -18.0)
                                    if event
                                        run_len[i, j] += 1
                                    else
                                        if run_len[i, j] >= min_consecutive_days
                                            cold_days_count[i, j] += run_len[i, j]
                                        end
                                        run_len[i, j] = 0
                                    end
                                end
                            end
                        end
                    end
                end
            end

            print_progress("Grand froid stream", step, total_steps)
        end
    end

    finalize_open_runs!(run_len, cold_days_count, min_consecutive_days)

    cold_hours_grid = fill(NaN, n_lon, n_lat)
    @inbounds for j in 1:n_lat
        for i in 1:n_lon
            if land_mask[i, j]
                cold_hours_grid[i, j] = cold_days_count[i, j] * day_hours
            end
        end
    end
    return cold_hours_grid
end


# ================================================================================================
# Visualisation
# ================================================================================================

"""
    vizumap(data_2d, weights_file; outline=true, ...)

Affiche une heatmap spatiale du signal fourni (`data_2d`) avec :
- masque France,
- contour de la France superposé (si `outline=true`).
"""
function vizumap(
    data_2d,
    weights_file;
    outline::Bool=true,
    title::AbstractString="Heures de grand froid",
    palette=:inferno,
    linecolor=:black,
    linewidth::Float64=1.2
)
    ds = NCDataset(weights_file)
    weights = ds["weights_frac"][:, :]
    lats = ds["latitude"][:]
    lons = ds["longitude"][:]
    close(ds)

    current_map = (ndims(data_2d) == 3) ? data_2d[:, :, 1]' : data_2d'
    mask = country_mask_for_map(weights, size(current_map))
    current_map[.!mask] .= NaN

    x_plot, y_plot, lon_idx, lat_idx = prepare_heatmap_axes(lons, lats, current_map)
    z_plot = current_map[lat_idx, lon_idx]

    valid_data = Float64[]
    for v in z_plot
        if !ismissing(v) && !isnan(Float64(v))
            push!(valid_data, Float64(v))
        end
    end

    isempty(valid_data) && error("La carte ne contient aucune valeur valide à afficher.")
    min_val, max_val = extrema(valid_data)

    p = heatmap(
        x_plot,
        y_plot,
        z_plot,
        title=title,
        clims=(min_val, max_val),
        c=palette,
        xlabel="Longitude",
        ylabel="Latitude",
        aspect_ratio=:equal,
        right_margin=5Plots.mm,
        yflip=false
    )

    if outline
        outline_segments = load_outline_segments(FRANCE_SHP)
        add_country_outline!(p, outline_segments; linecolor=linecolor, linewidth=linewidth)
    end

    return p
end


# ================================================================================================
# Sorties annuelles (GIF + courbe)
# ================================================================================================

"""
    compute_annual_cold_days_maps_streaming(...)

Calcule, pour chaque année, une carte du nombre de jours de grand froid par pixel.
Retourne :
- `annual_days_maps` de taille `[lon, lat, year]`
- `annual_mean_days` : moyenne des jours par pixel pour chaque année
"""
function compute_annual_cold_days_maps_streaming(
    data_folder::String,
    weights_file::String,
    year_range;
    cold_mode::Symbol=:orange,
    selected_months=collect(1:12),
    variable_name::String="t2m",
    min_consecutive_days::Int=1
)
    mode = normalize_cold_mode(cold_mode)
    years = year_range isa Integer ? [year_range] : collect(year_range)
    n_years = length(years)
    n_years == 0 && error("year_range est vide.")

    annual_days_maps = nothing
    annual_mean_days = Vector{Float64}(undef, n_years)

    for (k, y) in enumerate(years)
        yearly_hours = calculate_cold_days_streaming(
            data_folder,
            weights_file,
            y:y;
            cold_mode=mode,
            selected_months=selected_months,
            variable_name=variable_name,
            min_consecutive_days=min_consecutive_days
        )
        yearly_days = yearly_hours ./ 24.0

        if isnothing(annual_days_maps)
            annual_days_maps = Array{Float64}(undef, size(yearly_days, 1), size(yearly_days, 2), n_years)
        end
        annual_days_maps[:, :, k] = yearly_days

        total = 0.0
        n_valid = 0
        @inbounds for v in yearly_days
            if !isnan(v)
                total += v
                n_valid += 1
            end
        end
        annual_mean_days[k] = n_valid == 0 ? NaN : (total / n_valid)
        print_progress("Grand froid annuel", k, n_years)
    end

    return annual_days_maps, annual_mean_days, years
end

"""
    save_annual_maps_gif(...)

Crée un GIF (1 frame/an) à partir d'une série de cartes annuelles.
"""
function save_annual_maps_gif(
    annual_maps::Array{Float64, 3},
    years::AbstractVector{<:Integer},
    weights_file::String,
    gif_file::String;
    title_prefix::AbstractString="Grand froid (jours/pixel)",
    palette=:inferno,
    outline::Bool=true,
    fps::Int=4
)
    ds = NCDataset(weights_file)
    weights = ds["weights_frac"][:, :]
    lats = ds["latitude"][:]
    lons = ds["longitude"][:]
    close(ds)

    outline_segments = outline ? load_outline_segments(FRANCE_SHP) : Tuple{Vector{Float64}, Vector{Float64}}[]
    first_map = annual_maps[:, :, 1]'
    mask = country_mask_for_map(weights, size(first_map))
    x_plot, y_plot, lon_idx, lat_idx = prepare_heatmap_axes(lons, lats, first_map)

    valid_data = annual_maps[.!isnan.(annual_maps)]
    isempty(valid_data) && error("Aucune valeur valide dans annual_maps pour générer le GIF.")
    min_val, max_val = extrema(valid_data)

    anim = @animate for i in eachindex(years)
        current_map = annual_maps[:, :, i]'
        current_map[.!mask] .= NaN
        z_plot = current_map[lat_idx, lon_idx]

        p = heatmap(
            x_plot,
            y_plot,
            z_plot,
            title="$(title_prefix) - $(years[i])",
            clims=(min_val, max_val),
            c=palette,
            xlabel="Longitude",
            ylabel="Latitude",
            aspect_ratio=:equal,
            right_margin=5Plots.mm,
            yflip=false
        )
        if outline
            add_country_outline!(p, outline_segments)
        end
        p
    end

    gif(anim, gif_file, fps=fps)
    return gif_file
end

"""
    save_annual_curve(...)

Sauvegarde une courbe annuelle (moyenne des jours de grand froid par pixel).
"""
function save_annual_curve(
    years::AbstractVector{<:Integer},
    annual_mean_days::AbstractVector{<:Real},
    curve_file::String;
    title::AbstractString="Grand froid - jours moyens par pixel et par année",
    ylabel::AbstractString="Jours de grand froid / pixel"
)
    p = plot(
        collect(years),
        collect(annual_mean_days),
        linewidth=2.0,
        marker=:circle,
        markersize=3,
        xlabel="Année",
        ylabel=ylabel,
        title=title,
        legend=false
    )
    save_plot(p, curve_file)
    return p
end


# ================================================================================================
# Exécution script
# ================================================================================================

"""
    run_grand_froid_heatmap(; ...)

Pipeline complet prêt à l'emploi :
1. charge/agrège les températures horaires en Celsius,
2. calcule les heures de grand froid par pixel selon le mode choisi,
3. trace la heatmap avec contour France,
4. sauvegarde l'image.

Cette fonction est pratique quand le fichier est importé depuis un autre script.
"""
function run_grand_froid_heatmap(;
    data_folder::String=data_folder_temporel,
    weights_file::String=weight_temporel,
    year_range=1950:2025,
    selected_months=collect(1:12),
    output_file::Union{Nothing, String}=nothing,
    output_gif_file::Union{Nothing, String}=nothing,
    output_curve_file::Union{Nothing, String}=nothing,
    cold_mode::Symbol=:orange,
    outline::Bool=true,
    min_consecutive_days::Int=1,
    make_gif::Bool=true,
    make_curve::Bool=true
)
    mode = normalize_cold_mode(cold_mode)
    mode_txt = cold_mode_label(mode)
    output_file = isnothing(output_file) ? joinpath(plot_dir, "heatmap_grand_froid_1950_2025_$(mode_txt).png") : output_file

    println("Running heatmap pipeline...")
    println("Years: $(first(year_range)) -> $(last(year_range))")
    println("Mode grand froid: $mode_txt")
    annual_days_maps = nothing
    annual_mean_days = nothing
    years = nothing

    if make_gif || make_curve
        println("Step 1: Annual maps/curve data...")
        annual_days_maps, annual_mean_days, years = compute_annual_cold_days_maps_streaming(
            data_folder,
            weights_file,
            year_range;
            cold_mode=mode,
            selected_months=selected_months,
            variable_name="t2m",
            min_consecutive_days=min_consecutive_days
        )
        total_days_map = dropdims(sum(annual_days_maps, dims=3), dims=3)
        cold_hours = total_days_map .* 24.0
    else
        println("Step 1: Streaming grand-froid computation...")
        cold_hours = calculate_cold_days_streaming(
            data_folder,
            weights_file,
            year_range;
            cold_mode=mode,
            selected_months=selected_months,
            variable_name="t2m",
            min_consecutive_days=min_consecutive_days
        )
    end

    println("Step 2: Rendering map...")

    p = vizumap(
        cold_hours,
        weights_file;
        outline=outline,
        title="Heures de grand froid ($(first(year_range))-$(last(year_range))) - mode $(mode_txt)"
    )
    save_plot(p, output_file)
    println("Saved heatmap to: $output_file")

    if make_gif || make_curve
        if make_gif
            gif_file = isnothing(output_gif_file) ? joinpath(plot_dir, "heatmap_grand_froid_annual_$(mode_txt)_$(first(years))_$(last(years)).gif") : output_gif_file
            save_annual_maps_gif(
                annual_days_maps,
                years,
                weights_file,
                gif_file;
                title_prefix="Grand froid $(mode_txt) (jours/pixel)",
                outline=outline
            )
            println("Saved annual GIF to: $gif_file")
        end

        if make_curve
            curve_file = isnothing(output_curve_file) ? joinpath(plot_dir, "courbe_grand_froid_annual_$(mode_txt)_$(first(years))_$(last(years)).png") : output_curve_file
            save_annual_curve(
                years,
                annual_mean_days,
                curve_file;
                title="Grand froid $(mode_txt) - jours moyens par pixel et par année",
                ylabel="Jours de grand froid / pixel"
            )
            println("Saved annual curve to: $curve_file")
        end
    end

    return cold_hours, p
end


"""
    main(args=ARGS)

Point d'entrée en mode script.

Usage :
- `julia analyses_complementaires/France/heatmap_grand_froid_1950_2025.jl`
- `julia analyses_complementaires/France/heatmap_grand_froid_1950_2025.jl /chemin/sortie.png`
- `julia analyses_complementaires/France/heatmap_grand_froid_1950_2025.jl /chemin/sortie.png 1950 2025 rouge`
"""
function main(args=ARGS)
    output_file = length(args) >= 1 ? args[1] : nothing
    start_year = length(args) >= 2 ? parse(Int, args[2]) : 1950
    end_year = length(args) >= 3 ? parse(Int, args[3]) : 2025
    cold_mode = length(args) >= 4 ? Symbol(lowercase(args[4])) : :orange

    if end_year < start_year
        error("end_year ($end_year) doit être >= start_year ($start_year).")
    end

    run_grand_froid_heatmap(
        year_range=start_year:end_year,
        output_file=output_file,
        cold_mode=cold_mode
    )
    return nothing
end


# Équivalent Julia de `if __name__ == "__main__":`
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

#=
Mode script (défaut 1950-2025)
julia analyses_complementaires/France/heatmap_grand_froid_1950_2025.jl

Mode script avec sortie personnalisée
julia analyses_complementaires/France/heatmap_grand_froid_1950_2025.jl /tmp/grand_froid.png

Mode script avec intervalle d'années + mode (jaune|orange|rouge)
julia analyses_complementaires/France/heatmap_grand_froid_1950_2025.jl /tmp/grand_froid.png 1960 2025 rouge

Depuis un autre fichier :
include("analyses_complementaires/France/heatmap_grand_froid_1950_2025.jl")
cold_hours, p = run_grand_froid_heatmap(year_range=1950:2025, cold_mode=:orange)
=#
