using NCDatasets
using StatsPlots

const SAMPLE_NC = joinpath(@__DIR__, "..", "data", "raw-yearly-combined", "era5_land_hourly", "era5_land_fr_2025_12.nc")
const WEIGHTS_FOLDER = joinpath(@__DIR__, "..", "weights")
const PLOTS_FOLDER = joinpath(@__DIR__, "..", "plots")
# const FLIP_X = false
# const FLIP_Y = true
# const PLOT_TRANSPOSE = true
#
# function orient_map(a)
#     b = PLOT_TRANSPOSE ? permutedims(a, (2, 1)) : a
#     FLIP_Y && (b = reverse(b, dims=1))
#     FLIP_X && (b = reverse(b, dims=2))
#     return b
# end

function sanitize_filename(s::AbstractString)
    t = replace(s, r"[^\w\-]+" => "_")
    return strip(t, '_')
end

function align_temp_for_weights(temp, w)
    if size(temp) == size(w)
        return temp
    elseif size(temp) == reverse(size(w))
        return permutedims(temp, (2, 1))
    end
    error("Size mismatch: temp=$(size(temp)) weights=$(size(w))")
end

function find_coord_name(ds, keys)
    for k in keys
        if haskey(ds, k)
            return k
        end
    end
    error("coord not found among $(keys)")
end

function orient_for_plot(data, lons, lats)
    data_plot = data
    lons_plot = copy(lons)
    lats_plot = copy(lats)

    if length(lons_plot) == size(data_plot, 2) && length(lats_plot) == size(data_plot, 1)
        data_plot = permutedims(data_plot, (2, 1))
    end

    if lons_plot[1] > lons_plot[end]
        lons_plot = reverse(lons_plot)
        data_plot = reverse(data_plot, dims=1)
    end
    if lats_plot[1] > lats_plot[end]
        lats_plot = reverse(lats_plot)
        data_plot = reverse(data_plot, dims=2)
    end

    return data_plot, lons_plot, lats_plot
end

function plot_regions(; all::Bool=true, sample_nc::AbstractString=SAMPLE_NC,
                      weights_folder::AbstractString=WEIGHTS_FOLDER,
                      number::Int=1, use_cos::Bool=true)
    mkpath(PLOTS_FOLDER)

    ds_temp = NCDataset(sample_nc)
    lat_name = find_coord_name(ds_temp, ["latitude", "lat", "Latitude", "LAT"])
    lon_name = find_coord_name(ds_temp, ["longitude", "lon", "Longitude", "LON"])
    lats = ds_temp[lat_name][:]
    lons = ds_temp[lon_name][:]
    temp = ds_temp["t2m"][:, :, 1]
    temp = temp .- 273.15
    close(ds_temp)

    weight_var = use_cos ? "final_weights" : "weights_frac"

    if !all
        sample_wt = joinpath(weights_folder, "weights_region_$(number).nc")
        ds_w = NCDataset(sample_wt)
        w = ds_w[weight_var][:, :]
        region_name = get(ds_w.attrib, "region_name", string(number))
        close(ds_w)

        temp_aligned = align_temp_for_weights(temp, w)
        mask = w .< 0.5
        temp_masked = copy(temp_aligned)
        temp_masked[mask] .= NaN

        # heatmap(orient_map(temp_masked'), aspect_ratio=:equal, color=:thermal)
        plot_data, lons_plot, lats_plot = orient_for_plot(temp_masked, lons, lats)
        heatmap(lons_plot, lats_plot, permutedims(plot_data, (2, 1)),
                aspect_ratio=:equal, color=:thermal)
        region_slug = sanitize_filename(region_name)
        plot_file = joinpath(PLOTS_FOLDER, "plot_region_$(region_slug).png")
        isfile(plot_file) && rm(plot_file)
        savefig(plot_file)
        return plot_file
    end

    weight_files = sort(filter(f -> startswith(f, "weights_region_") && endswith(f, ".nc"),
                               readdir(weights_folder)))
    if isempty(weight_files)
        error("No weights_region_*.nc files found in $(weights_folder)")
    end

    for (idx, file) in enumerate(weight_files)
        ds_w = NCDataset(joinpath(weights_folder, file))
        w = ds_w[weight_var][:, :]
        region_name = get(ds_w.attrib, "region_name", string(idx))
        close(ds_w)

        temp_aligned = align_temp_for_weights(temp, w)
        mask = w .< 0.5
        temp_masked = copy(temp_aligned)
        temp_masked[mask] .= NaN

        # heatmap(orient_map(temp_masked'), aspect_ratio=:equal, color=:thermal)
        plot_data, lons_plot, lats_plot = orient_for_plot(temp_masked, lons, lats)
        heatmap(lons_plot, lats_plot, permutedims(plot_data, (2, 1)),
                aspect_ratio=:equal, color=:thermal)
        region_slug = sanitize_filename(region_name)
        plot_file = joinpath(PLOTS_FOLDER, "plot_region_$(region_slug).png")
        isfile(plot_file) && rm(plot_file)
        savefig(plot_file)
    end

    return nothing
end

plot_regions()
