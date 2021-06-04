#============================================================================
# nom : Fonction de découpage d'un raster
# auteur : Perle Charlot
# date de création : 2017-12-11
# dates de modification : 
#   2017- - :  
#----------------------------------------------------------------------------
# version de R : 2.15.2
# extension requises : 
library(sp)
library(rgdal)
library(raster)
#----------------------------------------------------------------------------
# Objectif : Découper un raster en suivant les contours de la rivière
# Entrées : Masque de découpage, raster à découper, satellite, année, mois
# Sorties : Raster (toutes les couches rasters sont rassemblées), avec seulement la rivière
#============================================================================


# **************
# * Constantes *
# **************

#chemin répertoire
setwd("C:/R_rep")

#importation du masque
masque <- readOGR("C:/R_rep/vectordata", "decoupe") 

#importation du raster
#bandes spectrales du domaine du visible : vert, bleu et rouge
green <- raster("rasterdata/decembre/pleiade_decembre2017_predecoupe.TIF", band = 1)
blue <- raster("rasterdata/decembre/pleiade_decembre2017_predecoupe.TIF", band = 2)
red <- raster("rasterdata/decembre/pleiade_decembre2017_predecoupe.TIF", band = 3)
#bande spectrale du domaine de l'infrarouge : proche infrarouge
pir <- raster("rasterdata/decembre/pleiade_decembre2017_pir_predecoupe.TIF")
  
# *************
# * Fonctions *
# *************

decoupageRaster <- function(masque, green, blue, red, pir){
  
  #Demande d'informations à l'utilisateur
  #pour pouvoir écrire le chemin d'enregistrement des fichiers
  satellite <- readline(prompt="Quel est le satellite utilisé ? ")
  annee <- readline(prompt="Quelle est l'année d'acquisition de la photographie ? ")
  mois <- readline(prompt="Quel est le mois d'acquisition de la photographie ? ")

  #vérification de l'égalité de la résolution
  #car bien souvent la résolution dans le PIR est meilleure que dans le visible
  res_visible <- res(green)
  res_pir <- res(pir)
  #si les résolutions sont différentes, on diminue la résolution la plus grande 
  if (res_pir[1] != res_visible[1]) {
    pir <- aggregate(pir, fact = 4, fun = mean)}
  
  #vérification que les 4 couches raster ont les mêmes dimensions
  row_green <- nrow(green)
  col_green <- ncol(green)
  row_pir <- nrow(pir)
  col_pir <- ncol(pir)
  
  #on vérifie que les 2 rasters aient le même nombre de lignes
  if (row_green != row_pir) {
    if(row_green < row_pir) {#si PIR a plus de lignes
      #transformation du raster PIR en matrice
      pir_ma <- as.matrix(pir)
      #suppression de la dernière ligne
      pir_ma <- pir_ma[-(row_pir),]
      #transformation matrice en raster
      vide <- raster(red)
      pir <- raster(pir_ma, xmn=0, xmx=col_green, ymn=0, ymx=row_green, template=vide)
      
       }
      #   else {#si GREEN a plus de lignes
      #   #transformation du raster GREEN en matrice
      #   green_ma <- as.matrix(green)
      #   #suppression de la dernière ligne
      #   green_ma <- pir_ma[-(row_green),]
      #   #transformation matrice en raster
      #   vide <- raster(red)
      #   green <- raster(green_ma, xmn=0, xmx=col_pir, ymn=0, ymx=row_pir, template=vide)}
    }
    
  #on vérifie que les 2 rasters aient le même nombre de colonnes
  if (col_green != col_pir) {
    if(col_green < col_pir) {#si PIR a plus de colonnes
      #transformation du raster PIR en matrice
      pir_ma <- as.matrix(pir)
      #suppression de la dernière ligne
      pir_ma <- pir_ma[,-(col_pir)]
      #transformation matrice en raster
      vide <- raster(red)
      pir <- raster(pir_ma, xmn=0, xmx=col_green, ymn=0, ymx=row_green, template=vide)
    }
      # else {#si GREEN a plus de colonnes
      # #transformation du raster GREEN en matrice
      # green_ma <- as.matrix(green)
      # #suppression de la dernière ligne
      # green_ma <- pir_ma[,-(col_green)]
      # #transformation matrice en raster
      # vide <- raster(red)
      # green <- raster(green_ma, xmn=0, xmx=col_pir, ymn=0, ymx=row_pir, template=vide)}
    }
  
  #vérification de l'égalité des SCR
  a <- identicalCRS(red, masque)
  #Si les SCR sont différents, on projète les coordonnées du vecteur masque dans le SCR du raster
  if ((a == "FALSE")) {masque<- spTransform(masque, CRS = crs(red))}
  
  #création des chemins d'enregistrement des rasters qui vont être découper
  chemin_green <- paste("rasterdata/", satellite, "_", mois, annee, "_green_decoupe.TIF", sep = "") 
  chemin_blue <- paste("rasterdata/", satellite, "_", mois, annee, "_blue_decoupe.TIF", sep = "") 
  chemin_red <- paste("rasterdata/", satellite, "_", mois, annee, "_red_decoupe.TIF", sep = "") 
  chemin_pir <- paste("rasterdata/", satellite, "_", mois, annee, "_pir_decoupe.TIF", sep = "") 
  
  #découpage des rasters + sauvegarde dans le répertoire
  green_decoupe <- mask(green, masque, filename = chemin_green, overwrite=TRUE)
  blue_decoupe <- mask(blue, masque, filename = chemin_blue, overwrite=TRUE)
  red_decoupe <- mask(red, masque, filename = chemin_red , overwrite=TRUE)
  pir_decoupe <- mask(pir, masque, filename = chemin_pir, overwrite=TRUE)
  
  #création NDVI = (PIR - RED) / (PIR + RED)
  ndvi <- (pir_decoupe - red_decoupe)/(pir_decoupe + red_decoupe)
  #réunion des 4 rasters découpé en un seul
  raster_total <- stack(green_decoupe, blue_decoupe, red_decoupe, pir_decoupe, ndvi)
  chemin_raster <- paste("rasterdata/", satellite, "_", mois, "_", annee, ".TIF", sep = "") 
  writeRaster(raster_total, chemin_raster, overwrite = TRUE)
  return(raster_total)
}
# ***********************
# * Programme principal *
# ***********************

raster <- decoupageRaster(masque, green, blue, red, pir)