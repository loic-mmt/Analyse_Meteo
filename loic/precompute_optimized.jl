using ArchGDAL
using Proj
using NCDatasets
using Statistics
using Base.Threads


# Build a Proj transformation from EPSG:4326 (lon/lat) to the shapefile CRS.
# This avoids relying on ArchGDAL's coordinate-transformation constructor, which changes across versions.
function proj_to_layer_crs(srs)
    crs_wkt = try
        ArchGDAL.toWKT(srs)
    catch
        ""
    end
    crs_str = !isempty(strip(crs_wkt)) ? crs_wkt : (try
        ArchGDAL.toPROJ4(srs)
    catch
        ""
    end)

    isempty(strip(crs_str)) && error("Layer CRS could not be exported to WKT/PROJ4. Is the .prj missing?")
    return Proj.Transformation("EPSG:4326", crs_str, always_xy=true)
end

# Iterate over all features of an OGR layer and return cloned geometries.
# ArchGDAL's feature-iteration API differs across versions, so we try multiple approaches.
function layer_feature_geoms(layer)
    geoms = ArchGDAL.IGeometry[]

    # 1) Newer ArchGDAL: ArchGDAL.eachfeature(layer) do feat ... end
    if isdefined(ArchGDAL, :eachfeature)
        eachf = getfield(ArchGDAL, :eachfeature)
        eachf(layer) do feat
            g = ArchGDAL.getgeom(feat)
            push!(geoms, ArchGDAL.clone(g))
        end
        return geoms
    end

    # 2) Some versions support iterating directly over the layer
    try
        for feat in layer
            g = ArchGDAL.getgeom(feat)
            push!(geoms, ArchGDAL.clone(g))
        end
        return geoms
    catch
        # fall through
    end

    # 3) Fallback: loop by feature index
    n = if isdefined(ArchGDAL, :nfeature)
        getfield(ArchGDAL, :nfeature)(layer)
    elseif isdefined(ArchGDAL, :getfeaturecount)
        getfield(ArchGDAL, :getfeaturecount)(layer)
    elseif isdefined(ArchGDAL, :featurecount)
        getfield(ArchGDAL, :featurecount)(layer)
    else
        error("Cannot count features on this ArchGDAL version (no eachfeature/iterator/featurecount).")
    end

    getfeat = if isdefined(ArchGDAL, :getfeature)
        getfield(ArchGDAL, :getfeature)
    elseif isdefined(ArchGDAL, :feature)
        getfield(ArchGDAL, :feature)
    else
        error("Cannot read features by index on this ArchGDAL version (no getfeature/feature).")
    end

    for idx in 0:(n - 1)
        feat = getfeat(layer, idx)
        g = ArchGDAL.getgeom(feat)
        push!(geoms, ArchGDAL.clone(g))
    end

    return geoms
end


function compute_bounds(coords::AbstractVector{<:Real})
    N = length(coords)
    N < 2 && error("Need >= 2 coords")
    mids = 0.5 .*(coords[1:end-1] .+ coords[2:end]) #  calcule les milieux entre chaque paire de coordonnées successives
    # coords[1:end-1] = tous sauf le dernier
    # coords[2:end] = tous sauf le premier
    b = Vector{Float64}(undef, N + 1) # N centres de pixels -> N+1 frontières.
    b[2:end-1] .= mids
    b[1] = coords[1] + (coords[1] - mids[1]) # premier bord
    b[end] = coords[end] + (coords[end] - mids[end]) # dernier bord
    return b
    
end

# ======== Inputs ========= 
# =========================

sample_nc = "../tests/era5_t2m_2025_01-02_fr.nc"     # un seul fichier pour la grille
shp = "../data/shapefiles/region.shp" 
out_weights_nc = "weights_france_final.nc"


# ======== grille NetCDF =========
# ================================

ds = NCDataset(sample_nc)

"""
Cherche le premier nom de variable présent dans 
le dataset parmi plusieurs possibilités -> pour la reproductibilité 
dans le cas ou les noms viendraient à changer.
"""
function find_coords(ds, keys)
    for k in keys
        if haskey(ds, k); return k; end
    end
    error("coord not found among $(keys)")
end

lat_name = find_coords(ds, ["latitude","lat","Latitude","LAT"])
lon_name = find_coords(ds, ["longitude","lon","Longitude","LON"])

# Lit toutes les valeurs de la variable latitude / longitude.
# NCDatasets peut retourner un Vector{Union{Missing,Float64}} (même si aucun missing n'est présent).
# -> On force ici des vecteurs Float64
lats_raw = ds[lat_name][:]
lons_raw = ds[lon_name][:]

lats = Float64.(coalesce.(lats_raw, NaN))
lons = Float64.(coalesce.(lons_raw, NaN))

if any(isnan, lats) || any(isnan, lons)
    error("Latitude/longitude coords contain missing values; cannot compute bounds safely.")
end

nlat = length(lats)
println("Grid size: nlat=", nlat, " nlon=", length(lons), " -> pixels=", nlat*length(lons))
nlon = length(lons)
println("Threads: ", Threads.nthreads())

# Calcule les bornes lat/lon.
lat_bnds = compute_bounds(lats)
lon_bnds = compute_bounds(lons)

close(ds)


# ======== Lecture du shapefile et union + reprojection WGS84 ======== 
# ====================================================================

fr_geom_src, pj_ll_to_src = ArchGDAL.read(shp) do dataset
    layer = ArchGDAL.getlayer(dataset, 0) # Prend la couche 0 (première couche du shapefile)

    # CRS du shapefile
    srs_src = ArchGDAL.getspatialref(layer)

    # Transformation (lon/lat EPSG:4326) -> CRS shapefile
    pj = proj_to_layer_crs(srs_src)

    geoms = layer_feature_geoms(layer)

    # Fusionne toutes les régions en un seul polygone union
    u = geoms[1]
    for k in 2:length(geoms)
        u = ArchGDAL.union(u, geoms[k])
    end

    return (u, pj)
end

println("Shapefile loaded and unioned. Starting weight computation...")

# Fast point-in-polygon test (lon/lat -> CRS shapefile -> contains)
@inline function inside_france_ll(fr_geom_src, pj_ll_to_src, lon::Float64, lat::Float64)::Bool
    X, Y = pj_ll_to_src(lon, lat)
    pt = ArchGDAL.createpoint(X, Y)
    return ArchGDAL.contains(fr_geom_src, pt)
end

# ======== fractions par pixel (subsampling) ===========
# ======================================================

k = 8       # 8x8 points par pixel (64 tests)-> 10 ou 12 pour précision
weights_frac = zeros(Float64, nlat, nlon) # matrice (nlat × nlon) initialisée à 0

#=
@threads sert à paralléliser une boucle for :
au lieu que le programme calcule i=1,2,3,... sur un seul cœur CPU,
Julia va répartir les itérations de la boucle sur plusieurs threads pour aller plus vite.
=#

# Progress tracking (thread-safe)
rows_done = Threads.Atomic{Int}(0)
step = max(1, Int(floor(0.05 * nlat)))  # print every ~5% of rows

@threads for i in 1:nlat

    # Per-row counters (avoid atomics in the hot inner loop)
    full0_row = 0
    full1_row = 0
    boundary_row = 0

    for j in 1:nlon
        # Bornes du pixel (i,j)
        lon0, lon1 = lon_bnds[j], lon_bnds[j+1]
        lat0, lat1 = lat_bnds[i], lat_bnds[i+1]

        lonm = 0.5 * (lon0 + lon1)
        latm = 0.5 * (lat0 + lat1)

        # 9 tests rapides (4 coins + 4 milieux d'arêtes + centre)
        any_in = false
        all_in = true

        # corners
        v = inside_france_ll(fr_geom_src, pj_ll_to_src, lon0, lat0); any_in |= v; all_in &= v
        v = inside_france_ll(fr_geom_src, pj_ll_to_src, lon1, lat0); any_in |= v; all_in &= v
        v = inside_france_ll(fr_geom_src, pj_ll_to_src, lon1, lat1); any_in |= v; all_in &= v
        v = inside_france_ll(fr_geom_src, pj_ll_to_src, lon0, lat1); any_in |= v; all_in &= v

        # edge midpoints
        v = inside_france_ll(fr_geom_src, pj_ll_to_src, lonm, lat0); any_in |= v; all_in &= v
        v = inside_france_ll(fr_geom_src, pj_ll_to_src, lonm, lat1); any_in |= v; all_in &= v
        v = inside_france_ll(fr_geom_src, pj_ll_to_src, lon0, latm); any_in |= v; all_in &= v
        v = inside_france_ll(fr_geom_src, pj_ll_to_src, lon1, latm); any_in |= v; all_in &= v

        # center
        v = inside_france_ll(fr_geom_src, pj_ll_to_src, lonm, latm); any_in |= v; all_in &= v

        if all_in
            # Pixel entièrement en France (très probable)
            weights_frac[i, j] = 1.0
            full1_row += 1
            continue
        elseif !any_in
            # Pixel entièrement hors France (très probable)
            weights_frac[i, j] = 0.0
            full0_row += 1
            continue
        end

        # Sinon: pixel de frontière -> subsampling k×k (inchangé)
        boundary_row += 1
        inside = 0
        total = k * k

        for a in 1:k
            x = lon0 + (a - 0.5) * (lon1 - lon0) / k
            for b in 1:k
                y = lat0 + (b - 0.5) * (lat1 - lat0) / k
                inside += inside_france_ll(fr_geom_src, pj_ll_to_src, x, y) ? 1 : 0
            end
        end

        weights_frac[i, j] = inside / total
    end

    # Optional per-row debug (comment out if too noisy)
    # println("Row ", i, ": full1=", full1_row, " full0=", full0_row, " boundary=", boundary_row)

    # update progress counter once per completed latitude row
    done = Threads.atomic_add!(rows_done, 1) + 1
    if done % step == 0 || done == nlat
        pct = round(100 * done / nlat; digits=1)
        println("Progress: ", done, "/", nlat, " rows (", pct, "%)")
    end
end

println("Weights computed. Saving to NetCDF...")

# ======== Poids finaux : fraction × cos(lat) =========
# =====================================================

"""
Convertit les latitudes en radians -> deg2rad.(lats) applique à chaque élément.
    • Calcule cos pour chaque latitude
    -> Parce qu’en lat/lon, la largeur réelle d’un degré de longitude diminue avec la latitude
"""
weights_lat = cos.(deg2rad.(lats)) # vecteur (nlat)
final_weights = weights_frac .* reshape(weights_lat, nlat, 1) # transforme en colonne (nlat×1) pour pouvoir multiplier chaque ligne de weights_frac

println("weights_frac min/max/mean = ",
        minimum(weights_frac), " / ", maximum(weights_frac), " / ", mean(weights_frac))
println("final_weights min/max/mean = ",
        minimum(final_weights), " / ", maximum(final_weights), " / ", mean(final_weights))


# ============ Sauvegarde dans un NetCDF ==============
# =====================================================

isfile(out_weights_nc) && rm(out_weights_nc)

wds = NCDataset(out_weights_nc, "c")
defDim(wds, lat_name, nlat)
defDim(wds, lon_name, nlon)

vlat = defVar(wds, lat_name, Float64, (lat_name,))
vlon = defVar(wds, lon_name, Float64, (lon_name,))
v1   = defVar(wds, "weights_frac",  Float64, (lat_name, lon_name))
v2   = defVar(wds, "final_weights", Float64, (lat_name, lon_name))

vlat[:] = lats # -> Écrit les données dans les variables.
vlon[:] = lons
v1[:, :] = weights_frac # [:, :] = toute la matrice.
v2[:, :] = final_weights

wds.attrib["description"] = "France weights: fraction-in-France (subsample) and final_weights=fraction*cos(lat)"
wds.attrib["subsample_k"] = string(k)

close(wds)

println("Saved: $out_weights_nc")

"""
Resultats > 25 min : 

weights_frac min/max/mean = 0.0 / 1.0 / 0.3868447580645161
final_weights min/max/mean = 0.0 / 0.7460573750616996 / 0.2658164279327024

Résultats Optimisés < 5.40 min : 
weights_frac min/max/mean = 0.0 / 1.0 / 0.38687406226556637
final_weights min/max/mean = 0.0 / 0.7460573750616996 / 0.26583931197628924

"""