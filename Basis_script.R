#BASIS SCRIPT_clean

#Prepare workspace ----
rm(big_object)       #remove large objects
gc()                 #free up the memory and report usage

#RUN FROM HERE
#install.packages(c("terra", "readr", "ape", "picante"))
library(terra)       #spatial data analysis
library(readr)       #read rectangular data (csv, tsv and fwf)
library(ape)         #read, manipulate and prune phylogenetic trees
library(picante)     #calculating PD
library(parallel)    #support for parallel computation, including by forking
library(foreach)     #support for parallel computation, execute repeated operations on multiple processors/cores
library(doParallel)  #a “parallel backend” for the foreach package, provides a mechanism needed to execute foreach loops in parallel
setwd("~/Desktop/Bachelor/Data bases") #sets working directory

#Global raster grid ----
# 1 degree resolution
tempwrlMap_1 <-rast(ncols=180*2, nrows=180) #creates a blank raster grid covering the globe.
#(ncols=180*2, nrows=180) → grid resolution (here: 360 x 180 = 1° resolution)
tempwrlMap_1[] <- 1 #fills all the cells with the value '1'
plot(tempwrlMap_1) #plots the empty grid

#PHYLACINE ----
##Raster stack all ----
#writeRaster(Phylacine_rast_all,"Phylacine_rast_all.tif") #Creates and saves the rasterstack as a tif file
Phylacine_rast_all <- rast("Phylacine_rast_all.tif") #Load tif file 

### 1° projection----
#Reprojects Phylacine raster to match the blank 1° world grid.
if(!"Phylacine_rast_all_proj.tif"%in%dir()){
  Phylacine_rast_all_proj <- project(Phylacine_rast_all, 
                                  tempwrlMap_1,
                                  #Projects raster to tempwrlMap (1°resolution)
                                  filename = "Phylacine_rast_all_proj.tif",
                                  #Saves projected raster as tif file
                                  overwrite = T)
                                  #Overrides excisting files with same name
} else {
  # Load the Reprojected Phylacine raster to match the blank 1° world grid.
  Phylacine_rast_all_proj <- rast("Phylacine_rast_all_proj.tif") #Load tif file
}
plot(sum(Phylacine_rast_all_proj), main ="All 1°") #adds up the raster → produces a species richness map (number of species present per cell), plots the result

### Clean up 'probability decimals' ----
#Gets rid of 'probability decimals' and creates a binary raster
  # 1 = species present
  # 0 = species absent 
if(!"Phylacine_rast_all_Poj_Binary.tif"%in%dir()){
  Phylacine_rast_all_Poj_Binary <- (Phylacine_rast_all_proj>0)*1
  #Convert all non-zero values (>0) to 1 to create a binary presence/absence raster stack
  plot(Phylacine_rast_all_Poj_Binary)
  #vizualize binary raster stack
  writeRaster(Phylacine_rast_all_Poj_Binary,
              "Phylacine_rast_all_Poj_Binary.tif",
              overwrite=T)
} else {
  Phylacine_rast_all_Poj_Binary <- rast("Phylacine_rast_all_Poj_Binary.tif") #Load tif file
}
plot(sum(Phylacine_rast_all_Poj_Binary), main ="All 1° binary") #adds up the raster → produces a species richness map (number of species present per cell), plots the result


##Trait data ----
#Trait_data <- read_csv("Phylacine 1.2.1/Data/Traits/Trait_data.csv") #import data
#saveRDS(Trait_data, "Trait_data.rds") # Save
Trait_data <- readRDS("Trait_data.rds") # Load back
#View(Trait_data) #Inspect data.frame
#names(Phylacine_rast_all_Poj_Binary) #View layer names in raster stack


###Filter Extinct species ----
#EP (extinct in prehistory, before 1500 CE)
#EX (extinct, after 1500 CE)
#EW (extinct in the wild)

#Remove extinct species
Trait_data_filtered <- Trait_data[!Trait_data$IUCN.Status.1.2 %in% c("EP", "EX", "EW"),]


###Habitat 'Vectors'----
#Overview of freshwater species (argument for classifying as terrestrial)
#TEST.FR <- Trait_data_filtered$Binomial.1.2[Trait_data_filtered$Freshwater == 1]
#Fresh.Trait_data.df <- subset(Trait_data_filtered, Binomial.1.2 %in% TEST.FR)
#unique(Fresh.Trait_data.df$Family.1.2) #30

#### Dual-habitat species ----
DHS.MA_TE <- Trait_data_filtered$Binomial.1.2[Trait_data_filtered$Marine == 1 & Trait_data_filtered$Terrestrial == 1] #Finds Dual-habitat species (Marine and Terrestrial)
DHS.MA_FR <- Trait_data_filtered$Binomial.1.2[Trait_data_filtered$Marine == 1 & Trait_data_filtered$Freshwater == 1] #Finds Dual-habitat species (Marine and Freshwater)
#DHS.MA_Ae <- Trait_data_filtered$Binomial.1.2[Trait_data_filtered$Marine == 1 & Trait_data_filtered$Aerial == 1] #Finds Dual-habitat species (0 found) (Marine and Aerial)
  #DHS.FR_TE_MA <- Trait_data_filtered$Binomial.1.2[Trait_data_filtered$Freshwater == 1 & Trait_data_filtered$Terrestrial == 1 & Trait_data_filtered$Marine == 1]
  #all(DHS.FR_TE_MA %in% DHS.MA_TE) #checks if all FR_TE_MA species are present in DHS.MA_TE (True)

Dual.habitat.species <- unique(c(DHS.MA_TE, DHS.MA_FR)) 

Dual.Trait_data <- subset(Trait_data_filtered, Binomial.1.2 %in% Dual.habitat.species) #Dataframe of Dual-habitat species (Marine and ...)
unique(Dual.Trait_data$Family.1.2) #Family overview of Dual-habitat species (11 different families)

#Selected families to be classified as Marine
Marine_families <- c("Otariidae","Delphinidae","Phocidae","Phocoenidae","Odobenidae","Trichechidae")

#Split Dual-habitat species into Marine vs Terrestrial based on family
Dual.Habitat_Marine <- Dual.Trait_data$Binomial.1.2[Dual.Trait_data$Family.1.2 %in% Marine_families]
  
Dual.Habitat_Terrestrial <- Dual.Trait_data$Binomial.1.2[!Dual.Trait_data$Family.1.2 %in% Marine_families]

##### Reclassification of Dual-habitat species----
#Reclassifies the chosen Dual-habitat species as only Marine (Otariidae, Delphinidae, Phocidae, Phocoenidae, Odobenidae, Trichechidae)
Trait_data_filtered[Trait_data_filtered$Binomial.1.2 %in% Dual.Habitat_Marine, c("Terrestrial", "Freshwater", "Aerial")] <-0  

#Reclassifies the chosen Dual-habitat species as not Marine (Terrestrial)(Mustelidae, Hippopotamidae, Muridae, Ursidae, Canidae)
Trait_data_filtered$Marine[Trait_data_filtered$Binomial.1.2 %in% Dual.Habitat_Terrestrial] <- 0 

#### Terrestrial Vector ----
#Terrestrial and/or Aerial and/or Freshwater
Terrestrial_species <- Trait_data_filtered$Binomial.1.2[Trait_data_filtered$Terrestrial == 1 | Trait_data_filtered$Aerial == 1 | Trait_data_filtered$Freshwater == 1] 

#### Marine Vector ----
Marine_species <- Trait_data_filtered$Binomial.1.2[Trait_data_filtered$Marine == 1]

intersect(Terrestrial_species, Marine_species)  #Safety check for double representation - (0)


###IUCN assessment 'Vectors' ----
#Creates habitat filtered lists only containing endangered species based on IUCN assessment
# CR (critically endangered)
# EN (endangered)
# VU (vulnerable)

#### Terrestrial Endangered Vector ----
above_EN_Terr <- Trait_data_filtered$Binomial.1.2[(Trait_data_filtered$IUCN.Status.1.2 %in% c("EN","CR")) &
                                           (Trait_data_filtered$Terrestrial == 1 | Trait_data_filtered$Aerial == 1| Trait_data_filtered$Freshwater == 1)]
#above_VU_Terr <- Trait_data_filtered$Binomial.1.2[(Trait_data_filtered$IUCN.Status.1.2 %in% c("VU","EN","CR")) &
#                                           (Trait_data_filtered$Terrestrial == 1 | Trait_data_filtered$Aerial == 1| Trait_data_filtered$Freshwater == 1)]

#### Marine Endangered Vector ----
above_EN_Mar <- Trait_data_filtered$Binomial.1.2[(Trait_data_filtered$IUCN.Status.1.2 %in% c("EN","CR")) &
                                           (Trait_data_filtered$Marine == 1)]
#above_VU_Mar <- Trait_data_filtered$Binomial.1.2[(Trait_data_filtered$IUCN.Status.1.2 %in% c("VU","EN","CR")) &
#                                           (Trait_data_filtered$Marine == 1)]


##Phylogeny ----
#Load Nexus file
#Small_trees <- read.nexus("./Phylacine 1.2.1/Data/Phylogenies/Small_phylogeny.nex")
Complete_trees <- read.nexus("./Phylacine 1.2.1/Data/Phylogenies/Complete_phylogeny.nex")

#Basic file check of myltiphylo's
#Small_trees    #multiphylo
#Complete_trees #multiphylo


#Remove large objects ----
#Removes excess dataframes, vectors and rasterstacks not needed for further analysis
rm(tempwrlMap_1,Phylacine_rast_all,Phylacine_rast_all_proj,Dual.Trait_data, Trait_data, DHS.MA_FR, DHS.MA_TE, Dual.Habitat_Marine, Dual.Habitat_Terrestrial, Dual.habitat.species, Marine_families)

#TO HERE


#References R ----
citation("terra")
citation("ape")
citation("picante")
citation("parallel")
citation("foreach")
citation("doParallel")
