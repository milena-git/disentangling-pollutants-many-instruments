setwd("Z:/Meteo LMD")
rm(list=ls())
install.packages("ncdf4")
library(ncdf4)
setwd("Données 06 mars 2019")

toread<-dir()[grepl(".nc$",dir())]

## Fichier exemple ##
u<-nc_open(toread[1])
u$dim$presnivs

# Pa
u$dim$presnivs$vals
u$dim$presnivs$vals[9]
#98.10^3 Pascal