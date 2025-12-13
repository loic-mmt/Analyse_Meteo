using Base.Threads

println("Threads disponibles = ", Threads.nthreads())

# Boucle parallèle : note que l'ordre d'affichage n'est pas garanti
@threads for i in 1:10
    # threadid() te dit quel thread exécute cette itération
    println("i=", i, " exécuté par thread ", threadid())
end