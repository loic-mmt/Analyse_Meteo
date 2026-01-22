# Analyse Météo — ERA5 France

Ce dépôt regroupe des scripts **Julia** et **Python** permettant de :
- télécharger des données **ERA5** (Copernicus CDS),
- construire un **masque / poids géographiques** sur la France,
- agréger des **climatologies** (mensuelles ou annuelles),
- produire des **visualisations** et analyses de tendance.

---

## Prérequis

- **Julia** (avec l’environnement du projet : `src/Project.toml`, `src/Manifest.toml`)
- **Python 3** + `cdsapi` (téléchargement ERA5)
- Dépendances système pour `ArchGDAL` / `LibGEOS` / `PROJ` (selon l’OS)
- Compte **Copernicus CDS** + clé API configurée (`~/.cdsapirc`)

---

## Installation (Julia)

Depuis la racine du dépôt :

```bash
julia --project=src -e 'using Pkg; Pkg.instantiate()'
```

---

## Téléchargement ERA5 (Python)

Le script `src/download.py` télécharge des fichiers NetCDF mensuels :
- zone : France (incluant la Corse)
- variable : température à 2 m (`t2m`)

Paramètres principaux à ajuster dans `src/download.py` :
- `START_YEAR` / `END_YEAR`
- `OUT_DIR` (par défaut : `era5_fr_t2m`)

Exécution :

```bash
python src/download.py
```

Recommandation : stocker les fichiers dans :

```
data/raw-yearly-combined/era5_fr_t2m
```

---

## Masque et poids France (Julia)

Le script `src/weights.jl` génère une matrice de poids à partir d’un shapefile France.
Il utilise un fichier ERA5 comme référence pour récupérer la grille (lon/lat).

Paramètres clés :
- `sample_nc` : fichier ERA5 local de référence
- `shp` : shapefile (par défaut : `data/shapefiles/region.shp`)
- `out_weights_nc` : fichier NetCDF de sortie

Exécution :

```bash
julia --project=src -e 'include("src/weights.jl")'
```

Le NetCDF produit contient :
- `weights_frac` : fraction du pixel dans la France
- `final_weights` : `weights_frac * cos(lat)` (pondération aire)

Un masque booléen peut être dérivé si nécessaire (ex. seuil `> 0.9`).

---

## Agrégation climatologique (Julia)

Le pipeline principal est dans `src/aggregate.jl` :
- lecture incrémentale des NetCDF (évite de charger tout en RAM),
- application du masque / des poids,
- conversion Kelvin → Celsius,
- export NetCDF final.

Exemples :

### Climatologie mensuelle

```julia
include("src/aggregate.jl")

cube = compute_general_climatology(
    "data/raw-yearly-combined/era5_fr_t2m",
    weights_bool_basic, # ou weights_prop_basic
    1950:2025;
    mode=:monthly,
    export_path="data/processed-means/mean_months_basic.nc",
)
```

### Climatologie annuelle

```julia
include("src/aggregate.jl")

cube = compute_general_climatology(
    "data/raw-yearly-combined/era5_fr_t2m",
    weights_bool_basic,
    1950:2025;
    mode=:yearly,
    export_path="data/processed-means/mean_years_basic.nc",
)
```

Notes :
- `weights_bool_basic`, `weights_prop_basic`, etc. sont chargés au début de `src/aggregate.jl`.
- si les fichiers de poids n’existent pas : corriger les chemins ou régénérer avec `weights.jl`.

---

## Visualisation et séries temporelles

Le script `src/visualization.jl` charge les climatologies et construit des séries temporelles France :
- `means_vector_calculation` : moyenne nationale par pas de temps
- tendances via `GLM` / `StatsPlots`

Exécution :

```bash
julia --project=src -e 'include("src/visualization.jl")'
```

---

## Cartes, animations et tendances spatiales

`src/Original.jl` contient des fonctions de référence (plus exploratoires) :
- cartes annuelles,
- tendances par pixel,
- animation GIF,
- p-values via `GLM`.

---

## Option performance : format JLD2

`src/JLD2 (in progress).jl` propose une conversion NetCDF → JLD2 pour accélérer les traitements.
Les chemins sont actuellement à adapter.

Exemple :

```julia
include("src/JLD2 (in progress).jl")
convert_to_fast_format("data/raw-yearly-combined/era5_fr_t2m", "data/fast/era5.jld2")
```

---

## Structure attendue

```
data/
  raw-yearly-combined/era5_fr_t2m/   # ERA5 mensuel (NetCDF)
  masks/                             # poids/masques France
  processed-means/                   # climatologies agrégées
  shapefiles/region.shp              # frontières France
src/
  *.jl                               # scripts Julia
  download.py                        # téléchargement ERA5
```

---

## Dépannage

- Échecs `ArchGDAL` / `LibGEOS` : vérifier l’installation **GDAL/GEOS/PROJ** sur la machine.
- Chemins codés en dur (`/mnt/...`, `/home/...`) : ajuster en tête de script.
- Variable NetCDF différente de `t2m` : modifier `variable_name` dans les fonctions.