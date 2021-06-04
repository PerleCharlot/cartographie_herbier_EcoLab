#Perle Charlot
#20/11/2017
#extraction des comportements spectraux sur image Pléiade
#chemin du répertoire
setwd("C:/R_rep")
#installation et chargement des packages
library(raster)
library(sp)
install.packages('rgdal')
library(rgdal)
library(lazyeval)
library(glue)
library(RStoolbox)
library(ggplot2)

#importation du raster, bande par bande
raster_total <- brick("rasterdata/decoupe_pleiade.TIF")
raster_B1 <- raster("rasterdata/decoupe_pleiade.TIF", band = 1)
raster_B2 <- raster("rasterdata/decoupe_pleiade.TIF", band = 2)
raster_B3 <- raster("rasterdata/decoupe_pleiade.TIF", band = 3)
raster_B4 <- raster("rasterdata/decoupe_pleiade.TIF", band = 4)

#extraction des valeurs radiométriques des pixels par bande spectrale
B1 <- getValues(raster_B1)
B2 <- getValues(raster_B2)
B3 <- getValues(raster_B3)
B4 <- getValues(raster_B4)

#création nouveau fichier des valeurs radiométriques pour chaque bande
#sans les 0 des contours blancs (diminution du poids donc temps de calcul)
nv <- c()
indice <- c()
for (i in 2:length(B1))
  if (xor(B1[i]>0, B1[i]==0 & B1[i+1]>0 & B1[i-1]>0)) {
    nv[i]=B1[i]
    indice[i]=i  }
valB1=nv[!is.na(nv)]
indice_pixel_B1=indice[!is.na(indice)]

#idem pour B2
nv <- c()
indice <- c()
for (i in 2:length(B2))
  if (xor(B2[i]>0, B2[i]==0 & B2[i+1]>0 & B2[i-1]>0)) {
    nv[i]=B2[i]
    indice[i]=i  }
valB2=nv[!is.na(nv)]
indice_pixel_B2=indice[!is.na(indice)]

#idem pour B3
nv <- c()
indice <- c()
for (i in 6:length(B3))
  if (xor(B3[i]>0, B3[i]==0 & B3[i+5]>0 & B3[i-5]>0)) {
    nv[i]=B3[i]
    indice[i]=i  }
valB3=nv[!is.na(nv)]
valB3 <- as.integer(valB3)
indice_pixel_B3=indice[!is.na(indice)]

#on vérifie que toutes mes bandes possèdent le même nombre de valeurs
a <- length(valB1)-length(valB3)
b <- length(valB1)-length(valB2)
c <- length(valB2)-length(valB3)
if (a != 0 ) {print("Les bandes 1 et 3 n'ont pas le même nombre de pixels")} else {print("Les bandes 1 et 3 ont le même nombre de pixels")}
if (b != 0 ) {print("Les bandes 1 et 2 n'ont pas le même nombre de pixels")} else {print("Les bandes 1 et 2 ont le même nombre de pixels")}
if (c != 0 ) {print("Les bandes 2 et 3 n'ont pas le même nombre de pixels")} else {print("Les bandes 2 et 3 ont le même nombre de pixels")}

#problème, B3 a moins de valeurs : des valeurs ont été oubliées (des 0 non contours)
#ajout manuel des 0 manquants

#fonction qui découpe en listes puis ajoute les 0 et recombine toutes 
#les listes en une seule

#les pixels manquants sont ceux apparaissants sur B1 mais pas sur B3
#il y en a 49
liste_manquant <- setdiff(indice_pixel_B1,indice_pixel_B3)
#liste par ordre décroissant
N <- sort(liste_manquant, decreasing = TRUE)
L1 <- valB3
L2 <- list()
liste_2 <- c()
liste_1 <- c()
for (i in N) {
  indice_matrice <- which(indice_pixel_B1 == i )
  L2 <- L1[indice_matrice:length(L1)]
  L1 <- L1[1:indice_matrice-1]
#rajouter le 0 manquant à la fin de la liste L2
  L2[length(L2)+1]=0
#liste_1 sert unique pour la dernière liste
  liste_1 <- append(liste_1, list(L1))
#liste_2 sert à stocker les listes intermédiaires 
  liste_2 <- append(liste_2, list(L2), after = 0) }

fin <- c()
for (i in 1:length(liste_2)){
  partie <- liste_2[i]
  fin <- append(fin, partie)
  L2_concat <- unlist(fin)  }

#extraction de la liste de liste_1 manquante
manq <- liste_1[length(liste_1)]
manq[length(manq)+1]=0
aj <- unlist(manq)
valeurs_de_B3 <- append(L2_concat, aj, after=0)
#enlever le dernier 0 
valB3 <- valeurs_de_B3[-738335]

#création matrice des indices et des valeurs
matrice_valeurs <- cbind (indice_pixel_B1, valB1, valB2, valB3)