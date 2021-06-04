#============================================================================
# nom : Affichage comportements spectraux des classes de classification
# auteur : Perle Charlot
# date de création : 2017-11-28
# dates de modification : 2017-12-12   
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
library(cowplot)
#----------------------------------------------------------------------------
# Objectif :Extraire les comportements spectraux des classes (valeurs radiométriques 
#           par bande) afin de savoir si nos classes définies sont distinctes les unes
#           des autres et respectivement homogènes
# Entrées : 
# Sorties : 
#============================================================================

# **************
# * Constantes *
# **************

#chemin répertoire
setwd("C:/R_rep")

#importation du raster, bande par bande
green <- raster("rasterdata/pleiade_septembre2017_green_decoupe.TIF")
blue <- raster("rasterdata/pleiade_septembre2017_blue_decoupe.TIF")
red <- raster("rasterdata/pleiade_septembre2017_red_decoupe.TIF")
pir <- raster("rasterdata/pleiade_septembre2017_pir_decoupe.TIF")

# *************
# * Fonctions *
# *************

comportementSpectral <- function(green, blue, red, pir){
  

#création NDVI = (PIR - RED) / (PIR + RED)
ndvi <- (pir - red)/(pir + red)

#réunion des 4 rasters découpé en un seul
raster_total <- stack(green, blue, red, pir,ndvi)

#importation des polygones d'entrainement
train_sept2017 <- readOGR("C:/R_rep/vectordata", "nouvelle_partition_v")

#vérification que les SCR sont égaux
a <- identicalCRS(raster_total, train_sept2017)
if ((a == "FALSE")) 
{train_sept2017<- spTransform(train_sept2017, CRS = crs(raster_total))}

#création vecteur unique classe de recouvrement
uniqueClasses <- unique(train_sept2017$Classe_rec)

#création de points d'échantillons au hasard, pris dans les polygones d'entrainements
set.seed(25)
for (i in 1:length(uniqueClasses)) {
  class_data <- subset(train_sept2017, Classe_rec == uniqueClasses[i])
  classpts <- spsample(class_data, type = "random", n=100)
  classpts$class <- rep(uniqueClasses[i], length(classpts))
  if (i == 1){xy <- classpts} else {xy <- rbind(xy, classpts)} }

#extraction des valeurs spectrales de chaque point d'échantillon
trainvals <- extract(raster_total, y= xy, cellnumbers = TRUE)
trainvals <- data.frame(response = xy$class, trainvals)

#élimination des points d'échantillonage ayant des pixels se chevauchants
mul <- any(duplicated(trainvals$cells))
if (mul == 'TRUE') {trainvals <- trainvals[!duplicated(trainvals$cells), -2]}

#affichage des coordonnées (valeurs radiométriques) des pixels
#appartenant à chaque classes

#création d'une palette de couleurs
#où le % d'herbiers est représenté par un dégradé de vert
#l'eau en bleu, et les galets en gris
palette <- c("#EDF8E9","#BAE4B3", "#74C476", "#31A354", "#006D2C","#0099CC","#CCCCCC")

#simplification d'écriture
Classes <- factor(trainvals$response)

#création du plot
R_PIR <- qplot(x = trainvals$pleiade_septembre2017_red_decoupe, 
               y = trainvals$pleiade_septembre2017_pir_decoupe, 
               data=trainvals, colour= Classes) + geom_point() +
  scale_color_manual(values=palette2) + 
  stat_ellipse() + theme_dark() + xlab("Rouge") + ylab("PIR") +
  ggtitle("Comportements spectraux") +labs(fill = "Classes") 

#affichage du plot
R_PIR
save_plot("rouge_pir.jpg", R_PIR)

#extraire la légende
legend <- get_legend(R_PIR)

G_PIR <- qplot(trainvals$pleiade_septembre2017_green_decoupe, 
               trainvals$pleiade_septembre2017_pir_decoupe , 
               data=trainvals, colour= Classes) +
  scale_color_manual(values=palette) + 
  stat_ellipse() + theme_dark() + xlab("Vert") + ylab("PIR") +
  ggtitle("Comportements spectraux") + labs(fill = "Classes") +
  theme(legend.position="none")

NDVI_R <- qplot(trainvals$pleiade_septembre2017_red_decoupe,
         trainvals$layer,data=trainvals, colour= Classes) +
  scale_color_manual(values=palette) + 
  stat_ellipse() + theme_dark() + xlab("Rouge") + ylab("NDVI") +
  ggtitle("Comportements spectraux") +
  theme(legend.position="none")

NDVI_G <- qplot(trainvals$pleiade_septembre2017_green_decoupe,
          trainvals$layer,data=trainvals, colour= Classes) +
  scale_color_manual(values=palette) +
  stat_ellipse() + theme_dark() + xlab("Vert") + ylab("NDVI") +
  ggtitle("Comportements spectraux") + theme(legend.position="none")

NDVI_B <- qplot(trainvals$pleiade_septembre2017_blue_decoupe,
          trainvals$layer,data=trainvals, colour= Classes) +
  scale_color_manual(values=palette) +
  stat_ellipse() + theme_dark() + xlab("Bleu") + ylab("NDVI") +
  ggtitle("Comportements spectraux") + theme(legend.position="none")

B_PIR <- qplot(trainvals$pleiade_septembre2017_blue_decoupe, 
               trainvals$pleiade_septembre2017_pir_decoupe , 
               data=trainvals, colour= Classes) + 
  scale_color_manual(values=palette) + 
  stat_ellipse() + theme_dark() + xlab("Bleu") + ylab("PIR") +
  ggtitle("Comportements spectraux") + labs(fill = "Classes") +
theme(legend.position="none")

B_G <- qplot(trainvals$pleiade_septembre2017_blue_decoupe, 
             trainvals$pleiade_septembre2017_green_decoupe , 
             data=trainvals, colour= Classes) + 
  scale_color_manual(values=palette) + 
  stat_ellipse() + theme_dark() + xlab("Bleu") + ylab("Vert") +
  ggtitle("Comportements spectraux") +labs(fill = "Classes") +
  theme(legend.position="none")
R_G <- qplot(trainvals$pleiade_septembre2017_red_decoupe, 
             trainvals$pleiade_septembre2017_green_decoupe , 
             data=trainvals, colour= Classes) +
  scale_color_manual(values=palette)+ 
  stat_ellipse() + theme_dark() + xlab("Rouge") + ylab("Vert") +
  ggtitle("Comportements spectraux") + labs(fill = "Classes") +
  theme(legend.position="none")
B_R <- qplot(trainvals$pleiade_septembre2017_blue_decoupe, 
             trainvals$pleiade_septembre2017_red_decoupe , 
             data=trainvals, colour= Classes) + 
  scale_color_manual(values=palette)+ 
  stat_ellipse() + theme_dark() + xlab("Bleu") + ylab("Rouge") +
  ggtitle("Comportements spectraux") +labs(fill = "Classes") +
  theme(legend.position="none")

neuf_plot <- grid.arrange(NDVI_G, NDVI_R,NDVI_B, G_PIR, R_PIR, B_PIR,
             B_R, R_G, B_G,ncol = 3, nrow = 3)
trois_plot <- grid.arrange(NDVI_R, R_PIR,R_G, ncol = 3, nrow=1)
grid.arrange(legend, ncol = 1 , nrow = 1)
}

# ***********************
# * Programme principal *
# ***********************
