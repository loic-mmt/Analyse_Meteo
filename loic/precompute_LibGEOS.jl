using ArchGDAL
using Proj
using NCDatasets
using Statistics
using Base.Threads
using LibGEOS

"""
Résumé script
	1.	Lit la grille lat/lon depuis un NetCDF ERA5 (centres de pixels).
	2.	Calcule les bornes de chaque pixel (rectangle) à partir des centres.
	3.	Lit le shapefile des régions, fait une union → une géométrie France.
	4.	Convertit la géométrie France en LibGEOS (plus rapide) et la prépare.
	5.	Calcul en 2 passes :
	    •	Pass 1 : 1 seul test “point dans France” au centre de chaque pixel → masque intérieur/extérieur.
	    •	Pass 2 : détecte les pixels de frontière (où ça change).
	    •	pixels “clairement dedans/dehors” → poids 1 ou 0 direct
	    •	pixels frontière → échantillonnage k*k (ex 8*8) pour estimer la fraction.
	6.	Multiplie par cos(lat) pour approximer l’aire.
	7.	Sauve dans un NetCDF.
"""


"""
Crée une transformation depuis EPSG:4326 (lon/lat) vers le CRS du shapefile.
Évite de dépendre du constructeur de transformation de coordonnées d'ArchGDAL, 
qui change d'une version à l'autre.

- try/catch: essaie d’exporter le CRS en WKT.
    - si ça échoue -> renvoie un string vide "".

- Si WKT non vide -> on l’utilise
    - sinon -> on essaie PROJ4
        - si ça échoue encore -> renvoie un string vide "".

- Si toujours vide → erreur

- always_xy=true -> garantit que l’ordre est bien (x=lon, y=lat).

"""
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


# ========== LibGEOS helpers ========== 
# =====================================


"""
Construire une géométrie GEOS à partir de WKT (le nom du constructeur diffère selon les versions de LibGEOS)

- wkt::AbstractString : typage -> wkt doit être un string.
- isdefined(module, :name) : check si une fonction existe dans cette version.

But : fonctionner même si la fonction a changé de nom selon la version.

"""
function geos_from_wkt(wkt::AbstractString)
    if isdefined(LibGEOS, :readgeom)
        return LibGEOS.readgeom(wkt)
    elseif isdefined(LibGEOS, :geomFromWKT)
        return LibGEOS.geomFromWKT(wkt)
    elseif isdefined(LibGEOS, :fromWKT)
        return LibGEOS.fromWKT(wkt)
    else
        error("LibGEOS: no WKT reader found (readgeom/geomFromWKT/fromWKT).")
    end
end


"""
Préparer une géométrie pour la répétition des taches (accélération)

Prepare geometry = LibGEOS crée une structure interne accélérant contains(...) répété.
- isdefined(module, :name) : check si une fonction existe dans cette version.

Si aucune fonction disponible -> renvoie g brut (moins rapide).
"""
function geos_prepare(g)
    if isdefined(LibGEOS, :prepare)
        return LibGEOS.prepare(g)
    elseif isdefined(LibGEOS, :preparegeom)
        return LibGEOS.preparegeom(g)
    elseif isdefined(LibGEOS, :prepared)
        return LibGEOS.prepared(g)
    else
        return g
    end
end

"""
Créer un point GEOS

- isdefined(module, :name) : check si une fonction existe dans cette version.

Si aucune fonction disponible -> on créé un WKT "POINT (x y)".
"""
function geos_point(x::Float64, y::Float64)
    if isdefined(LibGEOS, :Point)
        return LibGEOS.Point(x, y)
    elseif isdefined(LibGEOS, :createpoint)
        return LibGEOS.createpoint(x, y)
    elseif isdefined(LibGEOS, :point)
        return LibGEOS.point(x, y)
    else
        return geos_from_wkt("POINT ($(x) $(y))")
    end
end


#=
Vérifie si le point est dans la géométrie

@inline : demande au compilateur d’inliner (gain perf)
::Bool : garantit que la fonction renvoie un booléen

- isdefined(module, :name) : check si une fonction existe dans cette version.

=#

@inline function geos_contains(poly, pt)::Bool
    if isdefined(LibGEOS, :contains)
        return LibGEOS.contains(poly, pt)
    elseif isdefined(LibGEOS, :containsproperly)
        return LibGEOS.containsproperly(poly, pt)
    else
        error("LibGEOS: contains predicate not found.")
    end
end


"""
Lire les features du shapefile (API ArchGDAL variable) :

Crée un tableau vide qui contiendra des géométries GDAL (IGeometry)

Stratégie 1 - eachfeature : 
- getfield(Module, :name) récupère la fonction.
- do feat ... end -> on passe une fonction anonyme à eachf

Stratégie 2 - itéreration direct

Stratégie 3 - boucle par index
- récupère n (nb features) via plusieurs fonctions possibles
- récupère getfeat via plusieurs fonctions possibles
"""
function layer_feature_geoms(layer)
    geoms = ArchGDAL.IGeometry[]

    # Stratégie 1:
    if isdefined(ArchGDAL, :eachfeature)
        eachf = getfield(ArchGDAL, :eachfeature)
        eachf(layer) do feat
            g = ArchGDAL.getgeom(feat)
            push!(geoms, ArchGDAL.clone(g))
        end
        return geoms
    end

    # Stratégie 2:
    try
        for feat in layer
            g = ArchGDAL.getgeom(feat)
            push!(geoms, ArchGDAL.clone(g))
        end
        return geoms
    catch
    end

    # Stratégie 3:
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


"""
Calcul des bornes pixel.


"""
function compute_bounds(coords::AbstractVector{<:Real})
    N = length(coords)
    N < 2 && error("Need >= 2 coords")
    # coords[1:end-1] = tous sauf le dernier
    # coords[2:end] = tous sauf le premier
    mids = 0.5 .* (coords[1:end-1] .+ coords[2:end]) # calcule les milieux entre chaque paire de coordonnées successives
    b = Vector{Float64}(undef, N + 1) # N centres de pixels -> N+1 frontières.
    b[2:end-1] .= mids
    b[1] = coords[1] + (coords[1] - mids[1])         # premier bord
    b[end] = coords[end] + (coords[end] - mids[end]) # dernier bord
    return b
end

# ======== Inputs =========

sample_nc = "../tests/era5_t2m_2025_01-02_fr.nc"     # un seul fichier pour la grille
shp = "../data/shapefiles/region.shp"
out_weights_nc = "weights_france_final.nc"

# ======== grille NetCDF =========

ds = NCDataset(sample_nc)

"""
Cherche le premier nom de variable présent dans 
le dataset parmi plusieurs possibilités -> pour la reproductibilité 
dans le cas ou les noms viendraient à changer.
"""
function find_coords(ds, keys)
    for k in keys
        if haskey(ds, k)
            return k
        end
    end
    error("coord not found among $(keys)")
end

lat_name = find_coords(ds, ["latitude", "lat", "Latitude", "LAT"])
lon_name = find_coords(ds, ["longitude", "lon", "Longitude", "LON"])

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
nlon = length(lons)
println("Grid size: nlat=", nlat, " nlon=", nlon, " -> pixels=", nlat * nlon)
println("Threads: ", Threads.nthreads())

# Calcule les bornes lat/lon.
lat_bnds = compute_bounds(lats)
lon_bnds = compute_bounds(lons)

close(ds)

# ======== Lecture du shapefile (union) + proj lon/lat->CRS shapefile ========

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

# ======== Conversion en GEOS + prepared geometry ========
fr_wkt = ArchGDAL.toWKT(fr_geom_src) # WKT car GDAL géométrie != GEOS géométrie
fr_geos = geos_from_wkt(fr_wkt)
fr_geos_prepared = geos_prepare(fr_geos) # rapidité
println("LibGEOS geometry ready (prepared=", !(fr_geos_prepared === fr_geos), ")")
# `===` teste l'identité (même objet en mémoire), pas juste l'égalité des valeurs.
# Ici: si `fr_geos_prepared === fr_geos` alors la géométrie n'a pas pu être "préparée" (fallback), sinon elle l'a été.
# Le `!` devant un Bool inverse vrai/faux.

# Test d'appartenance d'un point à un polygone (longitude/latitude -> fichier de formes CRS -> contains)
@inline function inside_france_ll(fr_geos_prepared, pj_ll_to_src, lon::Float64, lat::Float64)::Bool
    # On part d'un point (lon, lat) en EPSG:4326 (WGS84).
    # Le shapefile peut être dans un autre CRS, donc on reprojette d'abord.
    X, Y = pj_ll_to_src(lon, lat)  # transforme (lon/lat) -> (X/Y) dans le CRS du shapefile

    # On construit ensuite un point GEOS (plus rapide pour répéter contains).
    # Float64(...) force le type (utile pour performance et compatibilité des constructeurs).
    pt = geos_point(Float64(X), Float64(Y))

    # Test géométrique: vrai si le point est dans le polygone France.
    return geos_contains(fr_geos_prepared, pt)
end

# ======== fractions par pixel (subsampling) =========

k = 8
weights_frac = zeros(Float64, nlat, nlon)


# ======== Pass 1: center mask (cheap) ========
# =============================================
center_in = falses(nlat, nlon)

println("Pass 1/2: computing center mask (1 contains per pixel)...")
# Idée: faire un premier passage très peu coûteux (1 test par pixel au centre).
# Ça donne une carte "dedans/dehors" (center_in) qu'on utilisera ensuite pour repérer uniquement les pixels proches de la frontière.
rows_done = Threads.Atomic{Int}(0)
step = max(1, Int(floor(0.05 * nlat))) # progrès : toutes les 5% de lignes

@threads for i in 1:nlat
    # `@threads` parallélise la boucle sur i : plusieurs lignes (latitudes) sont traitées en même temps.
    for j in 1:nlon
        lon0, lon1 = lon_bnds[j], lon_bnds[j+1] # récupère les bornes
        lat0, lat1 = lat_bnds[i], lat_bnds[i+1]
        lonm = 0.5 * (lon0 + lon1) # calcule le centre (lonm, latm)
        latm = 0.5 * (lat0 + lat1)
        center_in[i, j] = inside_france_ll(fr_geos_prepared, pj_ll_to_src, lonm, latm) # Test au centre.
    end

    done = Threads.atomic_add!(rows_done, 1) + 1
    if done % step == 0 || done == nlat
        pct = round(100 * done / nlat; digits=1)
        println("Pass 1 progress: ", done, "/", nlat, " rows (", pct, "%)")
    end
end


# ======== Pass 2: boundary detection + weights ========
# ======================================================
println("Pass 2/2: detecting boundary pixels and computing weights...")

# Identifier les pixels de bordure : un pixel est une bordure si l’un de ses 4 voisins diffère en center_in.
boundary = falses(nlat, nlon)

@threads for i in 1:nlat
    for j in 1:nlon
        v = center_in[i, j]  # statut du pixel au centre (true = dedans, false = dehors)
        b = false            # deviendra true si on détecte une discontinuité avec un voisin

        # Si un des 4 voisins (haut/bas/gauche/droite) a un statut différent, alors on est proche de la frontière.
        # Remarque syntaxe Julia: `if condition; action; end` permet d'écrire l'action sur la même ligne.
        if i > 1    && center_in[i-1, j] != v; b = true; end
        if i < nlat && center_in[i+1, j] != v; b = true; end
        if j > 1    && center_in[i, j-1] != v; b = true; end
        if j < nlon && center_in[i, j+1] != v; b = true; end

        boundary[i, j] = b
    end
end

# Correctif : dilater la limite de 1 pixel (quartier de 4 pixels)
# Ceci réduit les erreurs de classification des pixels proches de la limite tout en conservant k inchangé.
# Concrètement: si un pixel est marqué "boundary", on marque aussi ses 4 voisins.
# Objectif: inclure une petite zone tampon autour de la frontière pour éviter des pixels mal classés en intérieur/extérieur.
boundary_dil = copy(boundary)
@threads for i in 1:nlat
    for j in 1:nlon
        if boundary[i, j]
            boundary_dil[i, j] = true
            if i > 1;    boundary_dil[i-1, j] = true; end
            if i < nlat; boundary_dil[i+1, j] = true; end
            if j > 1;    boundary_dil[i, j-1] = true; end
            if j < nlon; boundary_dil[i, j+1] = true; end
        end
    end
end
boundary = boundary_dil

# Poids de remplissage : pixels intérieurs -> 0/1 directement, bordure -> calcul coûteux k×k
weights_frac .= 0.0 # remplit toute la matrice par 0 sans la recréer.

rows_done2 = Threads.Atomic{Int}(0) # compteur thread-safe pour afficher la progression de la passe 2
@threads for i in 1:nlat # utilisaion multi-coeur pour maximiser la vitesse d'execution
    for j in 1:nlon
        if !boundary[i, j]
            weights_frac[i, j] = center_in[i, j] ? 1.0 : 0.0 # Si pas frontière : poids direct 0 ou 1
            continue
        end

        # Pixel frontière -> calcul plus coûteux: on estime la fraction en échantillonnant k×k points dans le pixel.
        lon0, lon1 = lon_bnds[j], lon_bnds[j+1]  # bornes en longitude
        lat0, lat1 = lat_bnds[i], lat_bnds[i+1]  # bornes en latitude

        inside = 0
        total = k * k

        # On crée une grille régulière de k×k points.
        # (a - 0.5)/k place le point au centre de chaque "sous-cellule" -> échantillonnage plus stable.
        for a in 1:k
            x = lon0 + (a - 0.5) * (lon1 - lon0) / k
            for b in 1:k
                y = lat0 + (b - 0.5) * (lat1 - lat0) / k

                # `condition ? 1 : 0` = opérateur ternaire: ajoute 1 si dedans, sinon 0.
                inside += inside_france_ll(fr_geos_prepared, pj_ll_to_src, x, y) ? 1 : 0
            end
        end

        # Fraction estimée du pixel couverte par la France.
        weights_frac[i, j] = inside / total
    end

    done = Threads.atomic_add!(rows_done2, 1) + 1 # atomic_add! : évite que plusieurs threads cassent le compteur
    if done % step == 0 || done == nlat
        pct = round(100 * done / nlat; digits=1)
        println("Pass 2 progress: ", done, "/", nlat, " rows (", pct, "%)")
    end
end

nb_boundary = count(boundary)
println("Boundary pixels (dilated): ", nb_boundary, " / ", nlat*nlon, " (", round(100*nb_boundary/(nlat*nlon); digits=2), "%)")

println("Weights computed. Saving to NetCDF...")

# ======== Poids finaux : fraction × cos(lat) =========
weights_lat = cos.(deg2rad.(lats))
final_weights = weights_frac .* reshape(weights_lat, nlat, 1)

println("weights_frac min/max/mean = ",
        minimum(weights_frac), " / ", maximum(weights_frac), " / ", mean(weights_frac))
println("final_weights min/max/mean = ",
        minimum(final_weights), " / ", maximum(final_weights), " / ", mean(final_weights))

# ============ Sauvegarde dans un NetCDF ==============

isfile(out_weights_nc) && rm(out_weights_nc) # si le fichier existe déjà, on le supprime avant de le recréer

wds = NCDataset(out_weights_nc, "c")
defDim(wds, lat_name, nlat)
defDim(wds, lon_name, nlon)

vlat = defVar(wds, lat_name, Float64, (lat_name,))
vlon = defVar(wds, lon_name, Float64, (lon_name,))
v1   = defVar(wds, "weights_frac",  Float64, (lat_name, lon_name))
v2   = defVar(wds, "final_weights", Float64, (lat_name, lon_name))

vlat[:] = lats
vlon[:] = lons
v1[:, :] = weights_frac    # `[:, :]` = toute la matrice (toutes les lignes et toutes les colonnes)
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

Résultats Optimisés v2 < 4.11 min :
weights_frac min/max/mean = 0.0 / 1.0 / 0.3866279069767442
final_weights min/max/mean = 0.0 / 0.7460573750616996 / 0.2656593883012403

Résultats Optimisés v2 frontière dilatée < 7.12 min :
weights_frac min/max/mean = 0.0 / 1.0 / 0.3868447580645161
final_weights min/max/mean = 0.0 / 0.7460573750616996 / 0.2658164279327024

Résultats LibGEOS < 52 secondes :
weights_frac min/max/mean = 0.0 / 1.0 / 0.3868447580645161
final_weights min/max/mean = 0.0 / 0.7460573750616996 / 0.2658164279327024
"""