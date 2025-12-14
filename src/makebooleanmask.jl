using ArchGDAL
using Proj
using NCDatasets
using Base.Threads

# ==========================================
# 1. HELPER FUNCTIONS (Keep these for robustness)
# ==========================================
function proj_to_layer_crs(srs)
    wkt = try ArchGDAL.toWKT(srs) catch; "" end
    if isempty(wkt); wkt = ArchGDAL.toPROJ4(srs); end
    return Proj.Transformation("EPSG:4326", wkt, always_xy=true)
end

function get_union_geometry(layer)
    # Extracts all polygons and merges them into one (Union)
    geoms = ArchGDAL.IGeometry[]
    for i in 0:(ArchGDAL.nfeature(layer)-1)
        ArchGDAL.getfeature(layer, i) do feat
            push!(geoms, ArchGDAL.clone(ArchGDAL.getgeom(feat)))
        end
    end
    geom_union = geoms[1]
    for k in 2:length(geoms); geom_union = ArchGDAL.union(geom_union, geoms[k]); end
    return geom_union
end

# ==========================================
# 2. CONFIGURATION
# ==========================================
sample_nc     = "/mnt/data/ProjetMeteo/Analyse_Meteo/tests/era5_t2m_2025_01-02_fr.nc"
shp_file      = "/mnt/data/ProjetMeteo/Analyse_Meteo/data/shapefiles/region.shp"
output_file   = "mask_france_boolean.nc"

# ==========================================
# 3. LOAD DATA
# ==========================================
# A. Load Grid Coordinates
ds = NCDataset(sample_nc)
lats = Float64.(ds["latitude"][:]) # Adjust name if needed (e.g. "lat")
lons = Float64.(ds["longitude"][:])
close(ds)

nlat, nlon = length(lats), length(lons)
println("Grid: $nlat x $nlon")

# B. Load and Prepare Shapefile
shape_geom, projector = ArchGDAL.read(shp_file) do dataset
    layer = ArchGDAL.getlayer(dataset, 0)
    geom  = get_union_geometry(layer)       # Merge all regions
    srs   = ArchGDAL.getspatialref(layer)   # Get CRS
    trans = proj_to_layer_crs(srs)          # Build Proj transformer
    return (geom, trans)
end

# ==========================================
# 4. COMPUTE MASK (Center Point Check)
# ==========================================
println("Computing boolean mask...")
mask = zeros(Int8, nlat, nlon) # Int8 is sufficient (0 or 1)

@threads for i in 1:nlat
    for j in 1:nlon
        # 1. Get pixel center
        lat, lon = lats[i], lons[j]
        
        # 2. Project to Shapefile CRS (X, Y)
        X, Y = projector(lon, lat)
        
        # 3. Check "Point in Polygon"
        # createpoint is cheap, contains is the heavy check
        pt = ArchGDAL.createpoint(X, Y)
        if ArchGDAL.contains(shape_geom, pt)
            mask[i, j] = 1
        end
    end
end

# ==========================================
# 5. SAVE RESULT
# ==========================================
rm(output_file, force=true) # Delete if exists

NCDataset(output_file, "c") do ds
    defDim(ds, "latitude", nlat)
    defDim(ds, "longitude", nlon)

    # Save Coordinates
    defVar(ds, "latitude", lats, ("latitude",))
    defVar(ds, "longitude", lons, ("longitude",))
    
    # Save Mask
    v_mask = defVar(ds, "mask", Float64, ("latitude", "longitude"))
    v_mask[:, :] = mask
    
    v_mask.attrib["description"] = "Boolean Mask: 1=Inside, 0=Outside"
end

println("Success! Saved to $output_file")