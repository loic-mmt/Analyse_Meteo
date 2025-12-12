using ArchGDAL
using Proj
using NCDAtasets
using Statistics
using Base.Threads


function compute_bounds(coords::AbstractVector{<:Real})
    N = length(coords)
    N < 2 && error("Need >= coords")
    mids = 0.5 .*(coords[1:end-1] .+ coords[2:end])
    b = Vector{Float64}(undef, N + 1)
    b[2:end-1] .= mids
    b[1] = coords[end] + (coords[end] - mids[end])
    return b
    
end

# <----- Inputs -----> 

sample_nc = "era5_t2m_2025_01-02_fr.nc"     # un seul fichier pour la grille
shp = "../data/shapefiles/region.shp" 
out_weights_nc = "weights_france_final.nc"


# <----- NetCDF grid -----> 

ds = NCDataset(sample_nc)

function find_coords(ds, keys)
    for k in keys
        if haskey(ds, k); return k; end
    end
    error("coord not found among $(keys)")
end

lat_name = find_coord(ds, ["latitude","lat","Latitude","LAT"])
lon_name = find_coord(ds, ["longitude","lon","Longitude","LON"])

lats = ds[lat_name][:]
lons = ds[lon_name][:]

nlat = length(lats)
nlon = length(lats)

lat_bnds = compute_bounds(lats)
lon_bnds = compute_bounds(lons)

close(ds)


# <----- ShapeFile + union WGS84 -----> 

fr_geom_ll = ArchGDAL.read(shp) do dataset
    layer = ArchGDAL.getlayer(dataset, 0)
    srs_src = ArchGDAL.getspatialref(layer)
    srs_ll = ArchGDAL.importEPSG(4326)
    tr = ArchGDAL.createcoordinateTransformation(srs_src, srs_ll)

    geoms = ArchGDAL. IGeometry[]
    ArchGDAL.eachfeature(layer) do feat
        g = ArchGDAL.getgeom(feat)
        g2 = ArchGDAL.clone(g)
        ArchGDAL.transform!(g2, tr)
        push!(geoms, g2)
    end

    u = geoms[1]
    for k in 2:length(geoms)
        u = ArchGDAL.union(u, geoms[k])
    end
    u
end


# <----- fraction weights (subsampling) -----> 

k = 8       # 8x8 points par pixel (64 tests)-> 10 ou 12 pour précision
weights_frac = zeros(Float64, nlat, nlon)

@threads for i in 1:nlat
    for j in 1:nlat
        lon0, lon1 = lon_bnds[j], lon_bnds[j+1]
        lat0, lat1 = lat_bnds[i], lat_bnds[i+1]

        inside = 0
        total = k*k

        # grille régulière de points dans le pixel
        for a in 1:k
            x = lon0 + (a - 0.5) * (lon1 - lon0) / k
            for b in 1:k
                y =  lat0 + (b - 0.5) * (lat1 - lat0) / k
                pt = ArchGDAL.createpoint(x, y)
                inside += ArchGDAL.contains(fr_geom_ll, pt) ? 1 : 0
            end
        end

        weights_frac[i, j] = inside / total
    end
end


# <----- final weight (fraction * cos(lat) -> aire physique) -----> 

weights_lat = cos.(deg2rad.(lats))
final_weights = weights_frac .* reshape(weights_lat, nlat, 1)

println("weights_frac min/max/mean = ",
        minimum(weights_frac), " / ", maximum(weights_frac), " / ", mean(weights_frac))
println("final_weights min/max/mean = ",
        minimum(final_weights), " / ", maximum(final_weights), " / ", mean(final_weights))


# <----- Save weights to NetCDF (cache) -----> 

isfile(out_weights_nc) && rm(out_weights_nc)

wds = NCDataset(out_weights_nc, "c")
defDim(wds, lat_name, nlat)
defDim(wds, lon_name, nlon)

vlat = defVar(wds, lat_name, Float64, (lat_name,))
vlon = defVar(wds, lon_name, Float64, (lon_name,))
v1   = defVar(wds, "weights_frac",  Float64, (lat_name, lon_name))
v2   = defVar(wds, "final_weights", Float64, (lat_name, lon_name))

vlat[:] = lats
vlon[:] = lons
v1[:, :] = weights_frac
v2[:, :] = final_weights

wds.attrib["description"] = "France weights: fraction-in-France (subsample) and final_weights=fraction*cos(lat)"
wds.attrib["subsample_k"] = string(k)

close(wds)

println("Saved: $out_weights_nc")