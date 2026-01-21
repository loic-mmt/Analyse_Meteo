using NCDatasets
using StatsPlots

const SAMPLE_NC = joinpath(@__DIR__, "..", "data", "raw-yearly-combined", "era5_land_hourly", "era5_land_fr_2025_12.nc")
const WEIGHTS_FOLDER = joinpath(@__DIR__, "..", "data", "weights")
const PLOTS_FOLDER = joinpath(@__DIR__, "..", "plots")

function mask_from_weights(temp, w)
    if size(temp) == size(w)
        return w .< 0.5
    elseif size(temp) == reverse(size(w))
        return permutedims(w, (2, 1)) .< 0.5
    end
    error("Size mismatch: temp=$(size(temp)) weights=$(size(w))")
end

function plot_regions(; all::Bool=true, sample_nc::AbstractString=SAMPLE_NC,
                      weights_folder::AbstractString=WEIGHTS_FOLDER,
                      number::Int=1, use_cos::Bool=true)
    mkpath(PLOTS_FOLDER)

    ds_temp = NCDataset(sample_nc)
    temp = ds_temp["t2m"][:, :, 1]
    temp = temp .- 273.15
    close(ds_temp)

    weight_var = use_cos ? "final_weights" : "weights_frac"

    if !all
        sample_wt = joinpath(weights_folder, "weights_region_$(number).nc")
        ds_w = NCDataset(sample_wt)
        w = ds_w[weight_var][:, :]
        close(ds_w)

        mask = mask_from_weights(temp, w)
        temp_masked = copy(temp)
        temp_masked[mask] .= NaN

        heatmap(temp_masked', aspect_ratio=:equal, color=:thermal)
        plot_file = joinpath(PLOTS_FOLDER, "plot_region_$(number).png")
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
        close(ds_w)

        mask = mask_from_weights(temp, w)
        temp_masked = copy(temp)
        temp_masked[mask] .= NaN

        heatmap(temp_masked', aspect_ratio=:equal, color=:thermal)
        plot_file = joinpath(PLOTS_FOLDER, "plot_region_$(idx).png")
        isfile(plot_file) && rm(plot_file)
        savefig(plot_file)
    end

    return nothing
end

plot_regions()
