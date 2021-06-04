#============================================================================
# nom : Régression en classification
# auteur : Perle Charlot
# date de création : 2018-01-18
# dates de modification : 2018-01-19
#----------------------------------------------------------------------------
# version de R : 2.15.2
# extension requises : 
library(sp)
library(raster)
library(rgdal)
library(lazyeval)
library(lattice)
library(RStoolbox)
library(ggplot2)
library(RColorBrewer)
library(caret)
library(randomForest)
library(e1071)
library(nnet)
library(pls)
library(gbm)
library(caTools)
library(MASS)

# #############
# # Fonctions #
# #############

#Retourne la valeur la plus fréquente d'un vecteur
ClassePixel <- function (x){
  character <- as.character(x)
  if (is.na(character[1] == TRUE)){
    return(NA)
  }
  else {
    histo <- hist(x,breaks = c(0:7), plot = FALSE)
    classe <- histo$breaks
    classe <- classe[-1]
    eff <- histo$counts
    vl <- rbind(classe, eff)
    classe_maj <- which.is.max(vl[2,])
    return(classe_maj)
  }}

#séparation du jeu de données total en 2 jeux (entrainement et validation)
createSeparation <- function(data_classe, p){
  separation <- sample.split(data_classe, SplitRatio = p)
  return(separation)}

createTrainingDataSet <- function(data, separation){
  training_data <- subset(data, separation == TRUE)
  return(training_data)}

createValidationDataSet <- function(data, separation){
  validation_data <- subset(data, separation == FALSE)
  return(validation_data)}

#échantillonnage aléatoire de pixels au sein des polygones
SamplePixel <-  function(class, ensemble){
  for (i in 1:length(class)) {
    class_data <- subset(ensemble, Classe_rec == class[i])
    classpts <- spsample(class_data, type = "random", n=100)
    classpts$class <- rep(uniqueClasses[i], length(classpts))
    if (i ==1){xy <- classpts} else {xy <- rbind(xy, classpts)}}
  return(xy)}

# ##############
# # Constantes #
# ##############
setwd("C:/R_rep")

#importation des rasters, bande par bande
green <- raster("rasterdata/septembre/pleiade_septembre2017_green_decoupe.TIF")
blue <- raster("rasterdata/septembre/pleiade_septembre2017_blue_decoupe.TIF")
red <- raster("rasterdata/septembre/pleiade_septembre2017_red_decoupe.TIF")
pir <- raster("rasterdata/septembre/pleiade_septembre2017_pir_decoupe.TIF")
ndvi <- (pir - red)/(pir + red)

raster_total3 <- stack(green, blue, red)
raster_total5 <- stack(green, blue, red, pir,ndvi)

#importation du vecteur contenant les 60 placettes
total_data7 <- readOGR("C:/R_rep/vectordata", "nouvelle_partition")
total_data6 <- readOGR("C:/R_rep/vectordata", "6_classes")
total_data5 <- readOGR("C:/R_rep/vectordata", "5_classes")
total_data4 <- readOGR("C:/R_rep/vectordata", "4_classes")

#vérification des projections du raster et du vecteur : ils doivent être identiques
a <- identicalCRS(raster_total3, total_data7)
if ( (a == "FALSE") ) {
  total_data7 <- spTransform(total_data7, CRS = crs(raster_total3))} 
a <- identicalCRS(raster_total3, total_data6)
if ( (a == "FALSE") ) {
  total_data6 <- spTransform(total_data6, CRS = crs(raster_total3))}
a <- identicalCRS(raster_total3, total_data5)
if ( (a == "FALSE") ) {
  total_data5 <- spTransform(total_data5, CRS = crs(raster_total3))}
a <- identicalCRS(raster_total3, total_data4)
if ( (a == "FALSE") ) {
  total_data4 <- spTransform(total_data4, CRS = crs(raster_total3))}

#vérifier que les polygones soient bien dessinés
verif <- rgeos::gIsValid(total_data4)
if (verif == 'FALSE') {print("Redessinez les polygones")} else {print("Les polygones sont bien dessinés")}

#création de la palette de couleurs
pal7 <- c("#0099CC","#EDF8E9","#BAE4B3", "#74C476", "#31A354", "#006D2C","#CCCCCC")
pal6 <- c("#0099CC","#BAE4B3", "#74C476", "#31A354", "#006D2C","#CCCCCC")
pal5 <- c("#0099CC","#EDF8E9", "#74C476", "#006D2C","#CCCCCC")
pal4 <- c("#0099CC","#BAE4B3", "#006D2C","#CCCCCC")

iterations <- 10
#jeu de données
septembre <- total_data7$Classe_rec
#vecteur unitaire des classes
uniqueClasses <- unique(total_data$Classe_rec)

# #############
# # Programme #
# #############

# ******
# * RF *
# ******

# 1 : SUPERCLASS
kappa <- c()
overall_accuracy <- c()
dimensions <- dim(raster_total3)
nb_colonnes <- dimensions[1]*dimensions[2]
val_carte <- matrix(nrow = iterations, ncol = nb_colonnes)

#itérations de la classification
for (i in 1:iterations){
  classif <- superClass(raster_total3, model = "rf",
                        trainData = total_data6, trainPartition = 0.6,
                        responseCol = "Classe_rec", 
                        nSamples = 100, kfold = 3,
                        ntree = 1000,
                        mode = "classification")
  print(i)
  stats <- classif$validation$performance$overall
  overall_accuracy <- append(overall_accuracy, stats[1])
  kappa <- append(kappa, stats[2])
  map <- values(classif$map)
  val_carte <- rbind(val_carte, map)
  val_carte <- val_carte[-1,]
}
#Moyenne et écart type du kappa
moyenne_kappa <- mean(kappa)
moyenne_kappa
et_kappa <- sd(kappa)
et_kappa
moyenne_OA <- mean(overall_accuracy)
moyenne_OA*100
et_OA <- sd(overall_accuracy)
et_OA*100
freq_max_val <- apply(val_carte, 2, ClassePixel)
freq_max_val <- as.matrix(freq_max_val)
carte_freq_max <- raster(raster_total3)
carte_freq_max[] <- freq_max_val
plot(carte_freq_max, col=pal6)
writeRaster(filename = "results/test_classe/RF3_6-classes_SC.tif", carte_freq_max) 

#2 : Model + Predict
kappa <- c()
overall_accuracy <- c()
dimensions <- dim(raster_total3)
nb_colonnes <- dimensions[1]*dimensions[2]
val_carte <- matrix(nrow = iterations, ncol = nb_colonnes)
comp <- 0

for (j in 1:iterations){
  comp <- comp + 1
  print("Nous en sommes au tour n°")
  print(comp)
  sep <- createSeparation(septembre, 0.6)
  training_data <- createTrainingDataSet(total_data6, sep)
  validation_data <-createValidationDataSet(total_data6, sep)
  xy <- SamplePixel(uniqueClasses, training_data)
  trainvals <- extract(raster_total3, y= xy, cellnumbers = TRUE)
  trainvals <- data.frame(response = xy$class, trainvals)
  multiple <- any(duplicated(trainvals$cells))
  if (multiple == 'TRUE') {trainvals <- trainvals[!duplicated(trainvals$cells), -2]}
  model_1 <- randomForest(as.factor(response) ~ ., data = trainvals,na.action=na.omit, ntree = 1000)
  classif <- predict(raster_total3, model_1, type = "response", overwrite= TRUE)
  xy_val <- SamplePixel(uniqueClasses, validation_data)
  pred <- extract(classif, xy_val, cellnumbers = TRUE)
  dup <- duplicated(pred)
  pred <- pred[!dup, "layer"]
  obs <- xy_val$class[!dup]
  valFactor <- uniqueClasses[pred]
  statistiques_classification <- confusionMatrix(obs, reference = valFactor)
  stats <- statistiques_classification$overall
  overall_accuracy <- append(overall_accuracy, stats[1])
  kappa <- append(kappa, stats[2])
  map <- values(classif)
  val_carte <- rbind(val_carte, map)
  val_carte <- val_carte[-1,]}
# STATS  
moyenne_kappa <- mean(kappa)
moyenne_kappa
et_kappa <- sd(kappa)
et_kappa
moyenne_OA <- mean(overall_accuracy)
moyenne_OA*100
et_OA <- sd(overall_accuracy)
et_OA*100
freq_max_val <- apply(val_carte, 2, ClassePixel)
freq_max_val <- as.matrix(freq_max_val)
carte_freq_max <- raster(raster_total3)
carte_freq_max[] <- freq_max_val
plot(carte_freq_max, col=pal6)
carte_moy_RF <- writeRaster(filename = "results/test_classe/RF3_6-classes_MP.tif", carte_freq_max)

# *******
# * SVM *
# *******

# 1 : SUPERCLASS
kappa <- c()
overall_accuracy <- c()
dimensions <- dim(raster_total3)
nb_colonnes <- dimensions[1]*dimensions[2]
val_carte <- matrix(nrow = iterations, ncol = nb_colonnes)
#8 0.5
#itérations de la classification
for (i in 1:iterations){
  classif <- superClass(raster_total3, model = "svmRadial",
                        trainData = total_data6, trainPartition = 0.6,
                        responseCol = "Classe_rec", 
                        nSamples = 100, kfold = 3,
                        cost = 16, gamma = 4 ,
                        mode = "classification")
  print(i)
  stats <- classif$validation$performance$overall
  overall_accuracy <- append(overall_accuracy, stats[1])
  kappa <- append(kappa, stats[2])
  map <- values(classif$map)
  val_carte <- rbind(val_carte, map)
  val_carte <- val_carte[-1,]
}
#Moyenne et écart type du kappa
moyenne_kappa <- mean(kappa)
moyenne_kappa
et_kappa <- sd(kappa)
et_kappa
moyenne_OA <- mean(overall_accuracy)
moyenne_OA*100
et_OA <- sd(overall_accuracy)
et_OA*100
freq_max_val <- apply(val_carte, 2, ClassePixel)
freq_max_val <- as.matrix(freq_max_val)
carte_freq_max <- raster(raster_total3)
carte_freq_max[] <- freq_max_val
plot(carte_freq_max, col=pal6)
carte_moy_RF <- writeRaster(filename = "results/test_classe/SVM3_6-classes_SC.tif", carte_freq_max) 

#2 : Model + Predict
kappa <- c()
overall_accuracy <- c()
dimensions <- dim(raster_total3)
nb_colonnes <- dimensions[1]*dimensions[2]
val_carte <- matrix(nrow = iterations, ncol = nb_colonnes)
comp <- 0

for (j in 1:iterations){
  comp <- comp + 1
  print("Nous en sommes au tour n°")
  print(comp)
  sep <- createSeparation(septembre, 0.6)
  training_data <- createTrainingDataSet(total_data6, sep)
  validation_data <-createValidationDataSet(total_data6, sep)
  xy <- SamplePixel(uniqueClasses, training_data)
  trainvals <- extract(raster_total3, y= xy, cellnumbers = TRUE)
  trainvals <- data.frame(response = xy$class, trainvals)
  multiple <- any(duplicated(trainvals$cells))
  if (multiple == 'TRUE') {trainvals <- trainvals[!duplicated(trainvals$cells), -2]}
  model_1 <- svm(response ~ ., data = trainvals, kernel = "radial",type = "C-classification" ,cost = 16, gamma = 4 )
  classif <- predict(raster_total3, model_1, type = "response", overwrite= TRUE)
  xy_val <- SamplePixel(uniqueClasses, validation_data)
  pred <- extract(classif, xy_val, cellnumbers = TRUE)
  dup <- duplicated(pred)
  pred <- pred[!dup, "layer"]
  obs <- xy_val$class[!dup]
  valFactor <- uniqueClasses[pred]
  statistiques_classification <- confusionMatrix(obs, reference = valFactor)
  stats <- statistiques_classification$overall
  overall_accuracy <- append(overall_accuracy, stats[1])
  kappa <- append(kappa, stats[2])
  map <- values(classif)
  val_carte <- rbind(val_carte, map)
  val_carte <- val_carte[-1,]}

# STATS  
moyenne_kappa <- mean(kappa)
moyenne_kappa
et_kappa <- sd(kappa)
et_kappa
moyenne_OA <- mean(overall_accuracy)
moyenne_OA*100
et_OA <- sd(overall_accuracy)
et_OA*100
freq_max_val <- apply(val_carte, 2, ClassePixel)
freq_max_val <- as.matrix(freq_max_val)
carte_freq_max <- raster(raster_total3)
carte_freq_max[] <- freq_max_val
plot(carte_freq_max, col=pal6)
writeRaster(filename = "results/test_classe/SVM3_6-classes_MP.tif", carte_freq_max)

# *******
# * MLC *
# *******
kappa <- c()
overall_accuracy <- c()
dimensions <- dim(raster_total3)
nb_colonnes <- dimensions[1]*dimensions[2]
val_carte <- matrix(nrow = iterations, ncol = nb_colonnes)

#itérations de la classification
for (i in 1:iterations){
  classif <- superClass(raster_total3, model = "mlc",
                        trainData = total_data6, trainPartition = 0.6,
                        responseCol = "Classe_rec", 
                        nSamples = 100, kfold = 3,
                        mode = "classification")
  print(i)
  stats <- classif$validation$performance$overall
  overall_accuracy <- append(overall_accuracy, stats[1])
  kappa <- append(kappa, stats[2])
  map <- values(classif$map)
  val_carte <- rbind(val_carte, map)
  val_carte <- val_carte[-1,]#supprimer la ligne initiale vide
}
#Moyenne et écart type du kappa
moyenne_kappa <- mean(kappa)
moyenne_kappa
et_kappa <- sd(kappa)
et_kappa
moyenne_OA <- mean(overall_accuracy)
moyenne_OA*100
et_OA <- sd(overall_accuracy)
et_OA*100
freq_max_val <- apply(val_carte, 2, ClassePixel)
freq_max_val <- as.matrix(freq_max_val)
carte_freq_max <- raster(raster_total3)
carte_freq_max[] <- freq_max_val
plot(carte_freq_max, col=pal6)
carte_moy_RF <- writeRaster(filename = "results/test_classe/MLC3_6-classes_SC.tif", carte_freq_max) 

# *******
# * kNN *
# *******
kappa <- c()
overall_accuracy <- c()
dimensions <- dim(raster_total3)
nb_colonnes <- dimensions[1]*dimensions[2]
val_carte <- matrix(nrow = iterations, ncol = nb_colonnes)

#itérations de la classification
for (i in 1:iterations){
  classif <- superClass(raster_total3, model = "knn",
                        trainData = total_data6, trainPartition = 0.6,
                        responseCol = "Classe_rec", 
                        nSamples = 100, kfold = 3,
                        mode = "classification")
  print(i)
  stats <- classif$validation$performance$overall
  overall_accuracy <- append(overall_accuracy, stats[1])
  kappa <- append(kappa, stats[2])
  map <- values(classif$map)
  val_carte <- rbind(val_carte, map)
  val_carte <- val_carte[-1,]
}
#Moyenne et écart type du kappa
moyenne_kappa <- mean(kappa)
moyenne_kappa
et_kappa <- sd(kappa)
et_kappa
moyenne_OA <- mean(overall_accuracy)
moyenne_OA*100
et_OA <- sd(overall_accuracy)
et_OA*100
freq_max_val <- apply(val_carte, 2, ClassePixel)
freq_max_val <- as.matrix(freq_max_val)
carte_freq_max <- raster(raster_total3)
carte_freq_max[] <- freq_max_val
plot(carte_freq_max, col=pal6)
writeRaster(filename = "results/test_classe/kNN3_6-classes_SC.tif", carte_freq_max) 

# *******
# * ANN *
# *******
kappa <- c()
overall_accuracy <- c()
dimensions <- dim(raster_total5)
nb_colonnes <- dimensions[1]*dimensions[2]
val_carte <- matrix(nrow = iterations, ncol = nb_colonnes)

#itérations de la classification
for (i in 1:iterations){
  classif <- superClass(raster_total5, model = "nnet",
                        trainData = total_data6, trainPartition = 0.6,
                        responseCol = "Classe_rec", 
                        nSamples = 100, kfold = 3,
                        mode = "classification")
  print(i)
  stats <- classif$validation$performance$overall
  overall_accuracy <- append(overall_accuracy, stats[1])
  kappa <- append(kappa, stats[2])
  map <- values(classif$map)
  val_carte <- rbind(val_carte, map)
  val_carte <- val_carte[-1,]
}
#Moyenne et écart type du kappa
moyenne_kappa <- mean(kappa)
moyenne_kappa
et_kappa <- sd(kappa)
et_kappa
moyenne_OA <- mean(overall_accuracy)
moyenne_OA*100
et_OA <- sd(overall_accuracy)
et_OA*100
freq_max_val <- apply(val_carte, 2, ClassePixel)
freq_max_val <- as.matrix(freq_max_val)
carte_freq_max <- raster(raster_total5)
carte_freq_max[] <- freq_max_val
plot(carte_freq_max, col=pal6)
writeRaster(filename = "results/test_classe/ANN5_6-classes_SC.tif", carte_freq_max) 