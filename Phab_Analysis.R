library(readr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(lubridate)
library(psych)
library(FactoMineR)
library(ade4)
library(amap)
library(vegan)
library(knitr)
library(reshape2)
library(MASS)
library(data.table)
library(factoextra)

setwd("C:/SASAP_Bio/PCA")

#load data
Phab <- read_csv("C:/SASAP_Bio/PCA/Phab_HUClevel.csv")

#remove charcter vars and redundant habitat vars
Phab$name <- NULL
Phab$region <- NULL
Phab$id_numeric<- NULL
Phab$n_mines <- NULL
Phab$area_mines <- NULL
Phab$density_mi <- NULL
Phab$hfp_std <- NULL
Phab$forest <- NULL
Phab$hfp09_std <- NULL
Phab$region_id <- NULL
Phab$forest <- NULL
Phab$hfp09_mn <- NULL
Phab$std_l_area <- NULL

View(Phab)

# need to scale and center data
Phab_z<-scale(Phab, scale=TRUE, center=TRUE)
Phab_Z2 <- data.frame(t(na.omit(t(Phab_z))))

#following methodology from: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

#compute PCA - all data
res.pca <- prcomp(Phab_Z2, scale = FALSE)

#scree plot of eigenvalues - variance explained by each PC
fviz_eig(res.pca)

#graph of "individuals" - individuals with similar profiles grouped together - By regions? species? rather than hucs
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

#graph of variables positively correlated point to same side of plot - negatively correlated opposite side
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


#biplot of individuals and variables
fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)

#eigenvalues
eig.val <- get_eigenvalue(res.pca)
eig.val

# Results for Variables
res.var <- get_pca_var(res.pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 
# Results for individuals
res.ind <- get_pca_ind(res.pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 
