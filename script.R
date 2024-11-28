# Chargement des bibliothèques nécessaires
rm(list=ls())  
graphics.off() 
cat("\014")


library(coala)  # pour créer des modèles de coalescence
library(ape)    # pour manipuler et visualiser les arbres phylogénétiques
library(glue)
library(ggplot2)


################################################################################
# N CONSTANT ###################################################################
################################################################################

sample_size = 10
loci_length = 1000

# Définition du modèle de coalescence
model <- coal_model(
  sample_size = sample_size,         # N0 : Taille de l'échantillon de 10 individus
  loci_number = 1000,       # Nombre de loci (locus) dans le modèle
  loci_length = loci_length,       # Longueur du locus en base pairs (pb)
  ploidy = 2                # Ploïdie (nombre de copies de chaque gène par individu, ici diploïde)
)

# Specification du taux de mutation rate = 4 x N0 X u 
# Définition du taux de mutation (global pour le locus, à ajuster selon la longueur du locus)
model <- model + feat_mutation(rate = 1)  # Taux de mutation global pour chaque locus

# Ajout des statistiques résumées : Tajima's D (Dtaj), diversité nucléotidique (pi) et nombre de sites polymorphes (Sn)
model <- model + 
  sumstat_tajimas_d(name = "Dtaj") +            # Tajima's D pour la distribution des polymorphismes
  sumstat_nucleotide_div(name = "pi") +         # Diversité nucléotidique moyenne
  sumstat_seg_sites(name = "Sn")                # Nombre de sites polymorphes (0 ou 1)

model <- model + sumstat_trees()

# Simulation du modèle
sim_N_cst <- simulate(model, seed = 12)

# Lecture et visualisation de l'arbre phylogénétique simulé
tree <- read.tree(text = sim_N_cst$trees[[1]])
plot(tree)                             # Visualisation de l'arbre phylogénétique

# Calcul des statistiques résumées pour pi et Dtaj
mean_pi_N_cst = mean(sim_N_cst$pi, na.rm = TRUE)   # Moyenne de la diversité nucléotidique (pi)
hist_pi_N_cst = hist(sim_N_cst$pi, main = glue("mean = {mean_pi_N_cst}"),xlab="pi", breaks = "scott",col="#ADD8E6") # Histogramme de la diversité nucléotidique
abline(v=mean_pi_N_cst,col = "blue",lwd=3)

mean_Dtaj_N_cst = mean(sim_N_cst$Dtaj, na.rm = TRUE) # Moyenne de Tajima's D
hist_Dtaj_N_cst = hist(sim_N_cst$Dtaj, main = glue("mean = {mean_Dtaj_N_cst}"),xlab="Dtaj", breaks = "scott",col="#ADD8E6") # Histogramme de Tajima's D
abline(v=mean_Dtaj_N_cst,col = "blue",lwd=3)

# Accès aux sites polymorphes dans Sn
# sim$Sn [locus][[juste pour appel objet]]
sim_N_cst$Sn[1][[1]]  # site polymorphes du locus 1
sim_N_cst$Sn[2][[1]]  # site polymorphes du locus 1

# Nombre de sites polymorphes : obtenir la dimension de chaque matrice de sites
dim(sim_N_cst$Sn[1][[1]])[[2]]  # Nombre de sites polymorphes (colonnes) dans le premier locus


# CODE AJOUT RENAUD
# => Initialisation du SFS global
global_sfs_N_cst <- rep(0, sample_size + 1)  # Classes de fréquences de 0 à n

# Parcours de tous les loci simulés
for (i in seq_along(sim_N_cst$Sn)) {
  locus <- sim_N_cst$Sn[[i]]
  
  # Vérification si le locus contient des sites polymorphes
  if (!is.null(locus) && length(locus) > 0) {
    # Conversion en matrice
    locus_matrix <- as.matrix(locus)
    
    # Calcul des fréquences des allèles dérivés
    allele_frequencies <- colSums(locus_matrix)
    
    # Comptage des fréquences et mise à jour du SFS global
    locus_sfs <- table(factor(allele_frequencies, levels = 0:sample_size))
    global_sfs_N_cst <- global_sfs_N_cst + as.numeric(locus_sfs)
  }
}

names(global_sfs_N_cst) <- 0:sample_size

# Visualisation du SFS global
SFS_N_cst = barplot(
  global_sfs_N_cst,
  main = "Global Site Frequency Spectrum (Unfolded)",
  xlab = "Number of individuals with derived allele",
  ylab = "Number of polymorphic sites",
  col = "skyblue",
  names.arg = 0:sample_size,  # Affiche explicitement les 20 fréquences
)


















################################################################################
# Bottleneck ###################################################################
################################################################################

sample_size = 10
loci_length = 1000

# Définition du modèle de coalescence
model <- coal_model(
  sample_size = sample_size,         # N0 : Taille de l'échantillon de 10 individus
  loci_number = 1000,       # Nombre de loci (locus) dans le modèle
  loci_length = loci_length,       # Longueur du locus en base pairs (pb)
  ploidy = 2                # Ploïdie (nombre de copies de chaque gène par individu, ici diploïde)
)

# Specification du taux de mutation rate = 4 x N0 X u 
# Définition du taux de mutation (global pour le locus, à ajuster selon la longueur du locus)
model <- model + feat_mutation(rate = 1)  # Taux de mutation global pour chaque locus

# Changement de taille de la population à différents moments dans le temps
# Exemple : taille passée N0, puis 0.1xN0, puis N0
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
sim_N_var <- simulate(model, seed = 12)

# Lecture et visualisation de l'arbre phylogénétique simulé
tree <- read.tree(text = sim_N_var$trees[[1]])
plot(tree)                             # Visualisation de l'arbre phylogénétique

# Calcul des statistiques résumées pour pi et Dtaj
mean_pi_N_var = mean(sim_N_var$pi, na.rm = TRUE)   # Moyenne de la diversité nucléotidique (pi)
hist_pi_N_var = hist(sim_N_var$pi, main = glue("mean = {mean_pi_N_var}"),xlab="pi", breaks = "scott",col="#FFCCCB") # Histogramme de la diversité nucléotidique
abline(v=mean_pi_N_var,col = "red",lwd=3)

mean_Dtaj_N_var = mean(sim_N_var$Dtaj, na.rm = TRUE) # Moyenne de Tajima's D
hist_Dtaj_N_var = hist(sim_N_var$Dtaj, main = glue("mean = {mean_Dtaj_N_var}"),xlab="Dtaj", breaks = "scott",col="#FFCCCB") # Histogramme de Tajima's D
abline(v=mean_Dtaj_N_var,col = "red",lwd=3)

# Accès aux sites polymorphes dans Sn
# sim$Sn [locus][[juste pour appel objet]]
sim_N_var$Sn[1][[1]]  # site polymorphes du locus 1
sim_N_var$Sn[2][[1]]  # site polymorphes du locus 1

# Nombre de sites polymorphes : obtenir la dimension de chaque matrice de sites
dim(sim_N_var$Sn[1][[1]])[[2]]  # Nombre de sites polymorphes (colonnes) dans le premier locus


# CODE AJOUT RENAUD
# => Initialisation du SFS global
global_sfs_N_var <- rep(0, sample_size + 1)  # Classes de fréquences de 0 à n

# Parcours de tous les loci simulés
for (i in seq_along(sim_N_var$Sn)) {
  locus <- sim_N_var$Sn[[i]]
  
  # Vérification si le locus contient des sites polymorphes
  if (!is.null(locus) && length(locus) > 0) {
    # Conversion en matrice
    locus_matrix <- as.matrix(locus)
    
    # Calcul des fréquences des allèles dérivés
    allele_frequencies <- colSums(locus_matrix)
    
    # Comptage des fréquences et mise à jour du SFS global
    locus_sfs <- table(factor(allele_frequencies, levels = 0:sample_size))
    global_sfs_N_var <- global_sfs_N_var + as.numeric(locus_sfs)
  }
}

names(global_sfs_N_var) <- 0:sample_size

# Visualisation du SFS global
barplot(
  global_sfs_N_var,
  main = "Global Site Frequency Spectrum (Unfolded)",
  xlab = "Number of individuals with derived allele",
  ylab = "Number of polymorphic sites",
  col = "#FFCCCB",
  names.arg = 0:sample_size,  # Affiche explicitement les 20 fréquences
)

#### Comparaison avec N cst : 

gg = ggplot() +
  geom_histogram(aes(x = sim_N_cst$pi), 
                 bins = 60,  # ou scott
                 fill = "blue", color = "black", alpha = 0.4) +  
  geom_vline(aes(xintercept = mean_pi_N_cst), col = "blue", linetype = "dashed",lwd=1) +  
  
  geom_histogram(aes(x = sim_N_var$pi), 
                 bins = 60, # ou scott
                 fill = "red", color = "black", alpha = 0.4) +  #
  geom_vline(aes(xintercept = mean_pi_N_var), col = "red", linetype = "dashed",lwd=1) +  
  
  labs(title = glue("mean N cst = {mean_pi_N_cst}\nmean N var = {mean_pi_N_var}"),
       x = "pi", y = "Frequency")
plot(gg)

gg = ggplot() +
  geom_histogram(aes(x = sim_N_cst$Dtaj), 
                 bins = 30,  # ou scott
                 fill = "blue", color = "black", alpha = 0.4) +  
  geom_vline(aes(xintercept = mean_Dtaj_N_cst), col = "blue", linetype = "dashed",lwd=1) +  
  
  geom_histogram(aes(x = sim_N_var$Dtaj), 
                 bins = 30, # ou scott
                 fill = "red", color = "black", alpha = 0.4) +  #
  geom_vline(aes(xintercept = mean_Dtaj_N_var), col = "red", linetype = "dashed",lwd=1) +  
  
  labs(title = glue("mean N cst = {mean_Dtaj_N_cst}\nmean N var = {mean_Dtaj_N_var}"),
       x = "Dtaj", y = "Frequency")
plot(gg)

freq = c(0:10)
df1 <- data.frame(freq, global_sfs_N_cst,global_sfs_N_var)

# Créer le graphique
gg = ggplot(data = df1) +
  geom_col(aes(x = freq, y = global_sfs_N_cst), fill = "blue", alpha = 0.4, width = 0.8) +
  geom_col(aes(x = freq, y = global_sfs_N_var), fill = "red", alpha = 0.4, width = 0.8) +
  labs(title = "Global Site Frequency Spectrum (Unfolded)",
       x = "Number of individuals with derived allele",
       y = "Number of polymorphic sites") 
# Afficher le graphique
print(gg)












################################################################################
# MIGRATION ####################################################################
################################################################################

rm(list=ls())  
cat("\014")

loci_length = 1000

# Exemple avec deux populations
# Taille d'échantillon différente pour chaque population
model <- coal_model(
  sample_size = c(10, 8),    # Taille d'échantillon différente pour chaque population
  loci_number = 1000,        # Nombre de loci
  loci_length = loci_length,        # Longueur du locus
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

# Simulation de base
sim <- simulate(model)
fst_values <- 1 - (sim$pi1 + sim$pi2) / (2 * sim$pitot)

length(is.na(sim$pi1))
length(is.na(sim$pi2))
length(is.na(sim$pitot))
length(is.na(fst_values))

# pour stocker results
results <- data.frame(
  migration_rate = integer(0),
  pi1 = numeric(0),
  pi2 = numeric(0),
  pitot = numeric(0),
  fst = numeric(0)
)

# Boucle pour simuler plusieurs taux de migration et calculer Fst
for (migration_rate in c(0, 1, 5, 10, 50)) {
  # Modèle avec migration à taux variable
  model <- coal_model(
    sample_size = c(100, 80),   # Taille d'échantillon pour chaque population
    loci_number = 100,          # Nombre de loci
    loci_length = 100,          # Longueur du locus
    ploidy = 2                  # Ploïdie
  ) + 
    feat_mutation(rate = 1) + 
    feat_pop_merge(time = 1, pop_source = 1, pop_target = 2) + 
    feat_migration(rate = migration_rate, symmetric = T, time = 1) + # taux de migration = échange de n individus respectivements 
    
    sumstat_nucleotide_div(name = "pi1", population = 1) + 
    sumstat_nucleotide_div(name = "pi2", population = 2) + 
    sumstat_nucleotide_div(name = "pitot", population = "all")
  
  sim <- simulate(model)
  
  # Calcul de Fst
  fst_values <- 1 - (sim$pi1 + sim$pi2) / (2 * sim$pitot)
  print(mean(fst_values))
  

  # Ajouter les résultats au data frame
  results <- rbind(results, data.frame(
    migration_rate = migration_rate,
    pi1 = sim$pi1,
    pi2 = sim$pi2,
    pitot = sim$pitot,
    fst = fst_values
  ))
}



# 1. Histogrammes de Fst pour chaque taux de migration
gg <- ggplot(results, aes(x = fst, fill = factor(migration_rate))) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 60) +
  scale_fill_manual(values = c("blue", "green", "red", "purple", "orange")) + 
  labs(title = "Distribution Fst in fonction of migration rate",
       x = "Fst",
       y = "Frequency",
       fill = "Migration rate") +
  theme(legend.position = "top")
plot(gg)

# 2. Relation entre la diversité totale (pitot) et Fst pour chaque taux de migration
gg <- ggplot(results, aes(x = pitot, y = fst, color = factor(migration_rate))) +
  geom_point() +
  scale_color_manual(values = c("blue", "green", "red", "purple", "orange")) +
  labs(title = "Relation between Total diversity and Fst",
       x = "Total diversity pi total",
       y = "Fst",
       color = "Migration rate") 
  theme(legend.position = "top")
plot(gg)

