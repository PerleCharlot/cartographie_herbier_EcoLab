#============================================================================
# nom : Classification supervisée d'une image par l'algorithme Random Forest
# auteur : Perle Charlot
# date de création : 2017-11-29
# dates de modification : 2017-12-14   
#----------------------------------------------------------------------------
# version de R : 2.15.2
# extension requises : 
library(sp)
library(rgdal)
library(raster)
library(lazyeval)
library(lattice)
library(glue)
library(RStoolbox)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(randomForest)
library(cowplot)
#----------------------------------------------------------------------------
# Objectif : appliquer un algorithme de classification supervisée (RF) sur une image 
# Entrées : polygones d'entrainement (vecteur), polygones de validation (vecteur)
#           image à classifier (raster)
# Sorties : image classifiée (raster), statistiques de la classification (OA et Kappa)
#============================================================================

# **************
# * Constantes *
# **************

#chemin répertoire
setwd("C:/R_rep")

#importation raster, bande par bande
green <- raster("rasterdata/pleiade_septembre2017_green_decoupe.TIF")
blue <- raster("rasterdata/pleiade_septembre2017_blue_decoupe.TIF")
red <- raster("rasterdata/pleiade_septembre2017_red_decoupe.TIF")
pir <- raster("rasterdata/pleiade_septembre2017_pir_decoupe.TIF")
ndvi <- (pir - red)/(pir + red)
raster_total<- stack(green, blue, red, pir, ndvi)

#importation des vecteur d'entrainement et de validation
training_data <- readOGR("C:/R_rep/vectordata", "6_dec")
validation_data <- readOGR("C:/R_rep/vectordata", "polygone_validation_7c")

# *************
# * Fonctions *
# *************

classifRF <- function(training_data, validation_data, raster_total){
  
  #vérification des projections du raster et du vecteur : ils doivent être identiques
  a <- identicalCRS(raster_total, training_data)
  if ( (a == "FALSE") ) {
    training_data<- spTransform(training_data, CRS = crs(raster_total))} 
  b <- identicalCRS(raster_total, validation_data)
  if ( (b == "FALSE") ) {
    validation_data<- spTransform(validation_data, CRS = crs(raster_total))} 
  
  #vérifier que les polygones soient bien dessinés
  verif <- rgeos::gIsValid(training_data)
  try(if(verif == 'FALSE') stop("Redessinez les polygones"))
  verif2 <- rgeos::gIsValid(validation_data)
  try(if(verif2 == 'FALSE') stop("Redessinez les polygones"))
  
  #vérifier que les polygones d'entrainement et de validation possèdent 
  #une colonne commune identique, contenant le numéro des classes 
  #"Classe_rec"
  #col_recouv <- training_data$Classe_rec
  #try (if (col_recouv == "NULL") stop("Le vecteur d'entrainement ne possède pas de colonne 'Classe_rec' "))
  #col_recouv2 <- validation_data$Classe_rec
  #try (if (col_recouv2 == "NULL") stop("Le vecteur de validation ne possède pas de colonne 'Classe_rec' "))
  
  #lancement de la classification
  seilh_sept2017_sc <- superClass(img= raster_total, trainData= training_data,
                                  responseCol= "Classe_rec", model= "rf", mode="classification")
  
  #creation vecteur unique classe de recouvrement
  uniqueClasses <- unique(training_data$Classe_rec)
  
  #creation de points d'échantillons (pixels) au hasard, dans les polygones d'entrainement
  set.seed(25)
  for (i in 1:length(uniqueClasses)) {
    #choisir seulement les polygones des classes
    class_data <- subset(training_data, Classe_rec== uniqueClasses[i])
    #prendre des points aléatoires au sein de ces polygones
    classpts <- spsample(class_data, type = "random", n=100)
    #ajouter la colonne des classes dans un spatialspoints objet
    classpts$class <- rep(uniqueClasses[i], length(classpts))
    if (i ==1){xy <- classpts} else {xy <- rbind(xy, classpts)}
  }
  #extraction des valeurs spectrales de chaque point d'échantillon
  trainvals <- extract(raster_total, y= xy, cellnumbers = TRUE)
  trainvals <- data.frame(response = xy$class, trainvals)
  
  #regarder s'il y a des points en double
  multiple <- any(duplicated(trainvals$cells))
  # si oui, les supprimer
  if (multiple == 'TRUE') {trainvals <- trainvals[!duplicated(trainvals$cells), -2]}
  
  #lancer l'algorithme de classification RF
  #établissement de la relation entre les classes et leur reflectance selon chaque bande
  model_1 <- randomForest(response ~ ., data = trainvals,
                          na.action = na.omit, confusion = TRUE)
  
  #prédiction  du landcover grâce au model
  seilh_sept2017_sc <- predict(raster_total, model_1, 
                               filename= "results/seilh_sept2017_sc.TIF",format = "GTiff", 
                               datatype = "INT1U", type = "response", overwrite= TRUE)

  #count of pixel per classe
  seilh_sept2017_sc.freq <- freq(seilh_sept2017_sc, useNA = "no")
  
  #relié à la résolution spatial pour avoir surface
  resPleiade <- res(seilh_sept2017_sc)
  area_m2 <- seilh_sept2017_sc.freq[,"count"] * prod(resPleiade) 
  area <- data.frame(landcover = uniqueClasses, area_m2 = area_m2)
  area
  
  palette2 <- c("#0099CC","#EDF8E9","#BAE4B3", "#74C476", "#31A354", "#006D2C","#CCCCCC")
  carte <- ggR(seilh_sept2017_sc, geom_raster = TRUE)+
    scale_fill_gradientn(colours = palette2)
  carte
  return(seilh_sept2017_sc)

  #Evaluation précision de la classification
  #Classification avec le jeu de validation
  seilh_sept2017_sc <- superClass(raster_total, model = "rf",
                                  trainData = training_data,
                                  valData = validation_data,
                                  responseCol = "Classe_rec")
  seilh_sept2017_sc$validation$performance
  #Echantillonnage des points de validation
  set.seed(7)
  xy_val <- lapply(uniqueClasses, function(class_i){
    class_data <- subset(validation_data, Classe_rec == class_i)
    classpts <- spsample(class_data, type = "random", n= 100)
    classpts$class <- rep(class_i, length(classpts))
    return(classpts)
  })
  xy_val <- do.call("rbind", xy_val)
  pred <- extract(seilh_sept2017_sc$map, xy_val, cellnumbers = TRUE)
  dup <- duplicated(pred)
  pred <- pred[!dup, "Classe_rec"]
  obs <- xy_val$class[!dup]
  valFactor <- uniqueClasses[pred]
  
  #Afficher la matrice de confusion + coefficient Kappa de la classification
  confusionMatrix(obs, reference = valFactor)
  
  #Afficher la carte de classification
  palette2 <- c("#0099CC","#EDF8E9","#BAE4B3", "#74C476", "#31A354", "#006D2C","#CCCCCC")
  carte <- plot(seilh_sept2017_sc$map, col = palette2)
  seilh_sept2017_sc$model
  return(carte)
}

# ***********************
# * Programme principal *
# ***********************

classifRF(training_data, validation_data, raster_total)

####

library(party)
ctree_model <- ctree(response~., data = trainvals)
plot(ctree_model)
####