README - Analyse Meteo (ERA5 France)

But
Ce depot contient des scripts Julia et Python pour telecharger des donnees ERA5, calculer des masques geographiques sur la France, agregger des climatologies (mensuelles/annuelles), puis visualiser et analyser des tendances.

Guide rapide (tuto complet)

0) Prerequis
- Julia + packages du projet: `src/Project.toml` / `src/Manifest.toml`
- Python 3 avec `cdsapi` (pour le telechargement ERA5)
- Bibliotheques systeme pour ArchGDAL/LibGEOS/PROJ (selon votre OS)
- Acces Copernicus CDS + cle API configuree (~/.cdsapirc)

1) Installer l'environnement Julia
Depuis la racine du depot:
`julia --project=src -e 'using Pkg; Pkg.instantiate()'`

2) Telechargement des donnees ERA5 (Python)
Le script `src/download.py` telecharge des fichiers NetCDF mensuels:
- zone: France (incluant la Corse)
- variable: temperature 2m (t2m)

Dans `src/download.py`, ajuster:
- `START_YEAR` / `END_YEAR`
- `OUT_DIR` (par defaut `era5_fr_t2m`)

Execution:
`python src/download.py`

Recommandation:
- placer le resultat dans `data/raw-yearly-combined/era5_fr_t2m`
- ou modifier `OUT_DIR` pour ecrire directement dans ce dossier

3) Calcul des poids / masque France (Julia)
Le script `src/weights.jl` cree une matrice de poids a partir du shapefile France.
Il lit un fichier ERA5 pour recuperer la grille lon/lat.

Parametres importants (en haut du script):
- `sample_nc`: un fichier ERA5 local (ex: `data/raw-yearly-combined/era5_fr_t2m/era5_t2m_fr_1950_01.nc`)
- `shp`: shapefile (par defaut `data/shapefiles/region.shp`)
- `out_weights_nc`: fichier de sortie (ex: `data/masks/weights_prop_basic.nc`)

Execution:
`julia --project=src -e 'include("src/weights.jl")'`

Le fichier produit contient:
- `weights_frac`: fraction de pixel a l'interieur de la France
- `final_weights`: `weights_frac * cos(lat)` (pondere par l'aire)

Option: creer un masque booleen si besoin.
Exemple de regle utilisee dans les scripts: `weights .> 0.9` pour un masque 0/1.

4) Agregation climatologique (Julia)
Le pipeline principal est dans `src/aggregate.jl`:
- `accumulate_data` lit et accumule les NetCDF sans tout charger en RAM
- `finalize_cube` calcule les moyennes, applique le masque, convertit en Celsius
- `save_climatology_netcdf` exporte un NetCDF propre
- `compute_general_climatology` orchestre tout

Exemple mensuel (moyenne par mois):
`julia --project=src -q`
Puis dans le REPL Julia:
```
include("src/aggregate.jl")

cube = compute_general_climatology(
    "data/raw-yearly-combined/era5_fr_t2m",
    weights_bool_basic, # ou weights_prop_basic
    1950:2025;
    mode=:monthly,
    export_path="data/processed-means/mean_months_basic.nc",
)
```

Exemple annuel:
```
cube = compute_general_climatology(
    "data/raw-yearly-combined/era5_fr_t2m",
    weights_bool_basic,
    1950:2025;
    mode=:yearly,
    export_path="data/processed-means/mean_years_basic.nc",
)
```

Notes:
- `weights_bool_basic`, `weights_prop_basic`, etc. sont charges au debut de `src/aggregate.jl`
- Si ces fichiers n'existent pas, corriger les chemins ou regenerer les poids
- Les donnees sont supposees etre en Kelvin, conversion en Celsius faite dans le code

5) Visualisation / series temporelles
Le script `src/visualization.jl` charge les NetCDF agreges et calcule des series:
- `means_vector_calculation` calcule la temperature moyenne France par pas de temps
- utilisation de GLM/StatsPlots pour les tendances (voir `src/Original.jl`)

Execution (exemple):
`julia --project=src -e 'include("src/visualization.jl")'`

6) Cartes, animations et tendances spatiales
`src/Original.jl` contient des fonctions pour:
- cartes annuelles, tendances pixel, animation GIF
- calcul de p-values avec GLM
Ce fichier est plus experimental mais sert de reference pour reproduire les graphiques.

7) Option performance: format JLD2
`src/JLD2 (in progress).jl` propose une conversion des NetCDF vers un fichier JLD2 pour accelerer les calculs.
Les chemins sont actuellement hardcodes (a adapter). Utilisation typique:
```
include("src/JLD2 (in progress).jl")
convert_to_fast_format("data/raw-yearly-combined/era5_fr_t2m", "data/fast/era5.jld2")
```

Structure du dossier (attendue)
- `data/raw-yearly-combined/era5_fr_t2m/*.nc` -> fichiers ERA5 mensuels
- `data/masks/*.nc` -> masques/poids France
- `data/processed-means/*.nc` -> climatologies agreges
- `data/shapefiles/region.shp` -> frontieres France
- `src/*.jl` -> scripts Julia

Depannage
- Si `ArchGDAL` / `LibGEOS` echoue: verifier l'installation GDAL/GEOS/PROJ sur votre systeme.
- Si un script pointe vers `/mnt/...` ou `/home/...`: modifier les chemins en haut du fichier.
- Si la variable NetCDF n'est pas `t2m`: changer `variable_name` dans les fonctions.

Rappel: les scripts sont modulaires
Vous pouvez executer uniquement ce qui vous interesse:
- telecharger
- construire les poids
- agregger en NetCDF
- visualiser / analyser
