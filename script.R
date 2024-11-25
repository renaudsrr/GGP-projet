# Chargement des bibliothèques nécessaires
rm(list=ls())  
graphics.off() 
cat("\014")


library(coala)  # pour créer des modèles de coalescence
library(ape)    # pour manipuler et visualiser les arbres phylogénétiques
library(glue)



# Définition du modèle de coalescence
model <- coal_model(
  sample_size = 10,         # N0 : Taille de l'échantillon de 10 individus
  loci_number = 1000,       # Nombre de loci (locus) dans le modèle
  loci_length = 1000,       # Longueur du locus en base pairs (pb)
  ploidy = 2                # Ploïdie (nombre de copies de chaque gène par individu, ici diploïde)
)

# Définition du taux de mutation (global pour le locus, à ajuster selon la longueur du locus)
model <- model + feat_mutation(rate = 1)  # Taux de mutation global pour chaque locus

# Changement de taille de la population à différents moments dans le temps
# Exemple : taille passée 5xN0, puis 0.1xN0, et enfin à N0
model <- model +
  feat_size_change(new_size = 1, time = 1) +    # Taille de la population à 1xN0 à t=1
  feat_size_change(new_size = 0.1, time = 0.5)  # Taille de la population à 0.1xN0 à t=0.5

# Ajout des statistiques résumées : Tajima's D (Dtaj), diversité nucléotidique (pi) et nombre de sites polymorphes (Sn)
model <- model + 
  sumstat_tajimas_d(name = "Dtaj") +            # Tajima's D pour la distribution des polymorphismes
  sumstat_nucleotide_div(name = "pi") +         # Diversité nucléotidique moyenne
  sumstat_seg_sites(name = "Sn")                # Nombre de sites polymorphes (0 ou 1)

model <- model + sumstat_trees()

# Simulation du modèle
sim <- simulate(model)
                #, seed = 12)

# Lecture et visualisation de l'arbre phylogénétique simulé
tree <- read.tree(text = sim$trees[[1]])
plot(tree)                             # Visualisation de l'arbre phylogénétique

# Calcul des statistiques résumées pour pi et Dtaj
mean_pi = mean(sim$pi, na.rm = TRUE)   # Moyenne de la diversité nucléotidique (pi)
hist(sim$pi, main = glue("mean = {mean_pi}"),xlab="pi") # Histogramme de la diversité nucléotidique
abline(v=mean_pi,col = "red")

mean_Dtaj = mean(sim$Dtaj, na.rm = TRUE) # Moyenne de Tajima's D
hist(sim$Dtaj, main = glue("mean = {mean_Dtaj}"),xlab="Dtaj") # Histogramme de Tajima's D
abline(v=mean_Dtaj,col = "red")

# Accès aux sites polymorphes dans Sn
sim$Sn[1][[1]]  # Premier site polymorphe
sim$Sn[2][[1]]  # Deuxième site polymorphe

# Nombre de sites polymorphes : obtenir la dimension de chaque matrice de sites
dim(sim$Sn[1][[1]])[[2]]  # Nombre de sites polymorphes (colonnes) dans le premier locus

rm(list=ls())  
cat("\014")

# Exemple avec deux populations
# Taille d'échantillon différente pour chaque population
model <- coal_model(
  sample_size = c(10, 8),    # Taille d'échantillon différente pour chaque population
  loci_number = 1000,        # Nombre de loci
  loci_length = 1000,        # Longueur du locus
  ploidy = 2                 # Ploïdie (nombre de copies de chaque gène par individu)
)

model <- model + feat_mutation(rate = 1)

# Ajout d'événements de spéciation et de migration
model <- model + 
  feat_pop_merge(time = 1, pop_source = 1, pop_target = 2) +  # Fusion de populations à t=1
  feat_migration(rate = 1, symmetric = TRUE, time = 1)        # Migration symétrique entre populations à t=1

# Calcul de la diversité nucléotidique pour chaque population et pour la population totale
model <- model + 
  sumstat_nucleotide_div(name = "pi1", population = 1) +     # Diversité nucléotidique pour la population 1
  sumstat_nucleotide_div(name = "pi2", population = 2) +     # Diversité nucléotidique pour la population 2
  sumstat_nucleotide_div(name = "pitot", population = "all") # Diversité totale

# Simulation finale du modèle avec calcul du Fst
sim <- simulate(model, seed = 12)

# Boucle pour simuler plusieurs taux de migration et calculer Fst
for (migration_rate in c(0, 1, 5, 10, 50)) {
  # Modèle avec migration à taux variable
  model <- coal_model(
    sample_size = c(100, 80),   # Taille d'échantillon pour chaque population
    loci_number = 100,          # Nombre de loci
    loci_length = 100,          # Longueur du locus
    ploidy = 2                  # Ploïdie
  ) + feat_mutation(rate = 1) + 
    feat_pop_merge(time = 1, pop_source = 1, pop_target = 2) + 
    feat_migration(rate = migration_rate, symmetric = TRUE, time = 1) + 
    sumstat_nucleotide_div(name = "pi1", population = 1) + 
    sumstat_nucleotide_div(name = "pi2", population = 2) + 
    sumstat_nucleotide_div(name = "pitot", population = "all")
  
  # Calcul de Fst : Fst = 1 - (pi1 + pi2) / (2 * pitot)
  sim$fst <- 1 - (sim$pi1 + sim$pi2) / (2 * sim$pitot)
  
  # Simulation du modèle et affichage de Fst moyen
  sim <- simulate(model)
  print(mean(sim$fst))
}



# Calcul de Fst à partir de pi pour chaque population
sim$fst <- 1 - (sim$pi1 + sim$pi2) / (2 * sim$pitot)

# Visualisation du Fst et de la diversité totale
hist(sim$fst)                    # Histogramme de Fst
plot(sim$pitot, sim$fst)         # Relation entre la diversité totale et Fst
