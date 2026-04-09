# 1. Configuration des chemins et chargement
cd(joinpath(@__DIR__, "..", ".."))
include("../../demo_france/function2charge.jl")
include("../../demo_france/dataset2load.jl")

# Dossier de sortie spécifique
plot_dir = joinpath(@__DIR__, "plot")
isdir(plot_dir) || mkpath(plot_dir)

