#Data visualization 

#Prepare workspace ----
#install.packages(c("ggplot2","terra","tidyterra","rnaturalearth","rnaturalearthdata"))
library(ggplot2)            #Advanced data visualization and plotting
library(patchwork)          #Combine multiple ggplot panels into single figures
library(hexbin)             #Hexagonal binning for density plots and heatmaps
library(terra)              #Raster data manipulation and spatial analysis
library(tidyterra)          #Enables ggplot2 integration with terra rasters (geom_spatraster)
library(sf)                 #Simple features for vector spatial data (points, polygons, lines)
library(rnaturalearth)      #Access to Natural Earth vector map data
library(rnaturalearthdata)  #Provides country/continent boundary shapefiles for mapping

#World Map for ggplot ----
world <- ne_countries(scale = "medium", returnclass = "sf")

#Result df's----
##Terrestrial
if(!"Results_terr_df.rds" %in% dir()){
  results_terr <- read.csv("SES_Summ_Final_Terr.csv")
  saveRDS(results_terr, file = "Results_terr_df.rds") #saves df as rds file
} else {
  results_terr <- readRDS("Results_terr_df.rds") #Load rds file
}

##Marine
if(!"Results_mar_df.rds" %in% dir()){
  results_mar <- read.csv("SES_Summ_Final_Mar.csv")
  saveRDS(results_mar, file = "Results_mar_df.rds") #saves df as rds file
} else {
  results_mar <- readRDS("Results_mar_df.rds") #Load rds file
}

##Categorical comparison of results ----
#Create cell categories
results_terr$Category <- cut(results_terr$SES.Mean, 
                             breaks = c(-Inf, -2, 0, 2, Inf), #Define breaks ("boxes")
                             labels = c("Strongly Redundant", "Slightly Redundant", 
                                        "Slightly Unique", "Strongly Unique")) #Name breaks

results_mar$Category <- cut(results_mar$SES.Mean,
                            breaks = c(-Inf, -2, 0, 2, Inf), #Define breaks
                            labels = c("Strongly Redundant", "Slightly Redundant",
                                       "Slightly Unique", "Strongly Unique"))
#Count by category
#table(results_terr$Category)
#table(results_mar$Category)

###Cell distribution df (Table 1) ----
Cell_distribution <- data.frame(
  Category = c("Strongly Redundant", "Slightly Redundant", 
               "Slightly Unique", "Strongly Unique"),
  Terrestrial = c(
    sum(results_terr$Category == "Strongly Redundant"),
    sum(results_terr$Category == "Slightly Redundant"),
    sum(results_terr$Category == "Slightly Unique"),
    sum(results_terr$Category == "Strongly Unique")
  ),
  Proportion_Terr = c(((sum(results_terr$Category == "Strongly Redundant")/2514)*100), 
                 ((sum(results_terr$Category == "Slightly Redundant")/2514)*100),
                 ((sum(results_terr$Category == "Slightly Unique")/2514)*100),
                 ((sum(results_terr$Category == "Strongly Unique")/2514)*100)),
  Marine = c(
    sum(results_mar$Category == "Strongly Redundant"),
    sum(results_mar$Category == "Slightly Redundant"),
    sum(results_mar$Category == "Slightly Unique"),
    sum(results_mar$Category == "Strongly Unique")
  ),
  Proportion_Mar = c(((sum(results_mar$Category == "Strongly Redundant")/21819)*100), 
                      ((sum(results_mar$Category == "Slightly Redundant")/21819)*100),
                      ((sum(results_mar$Category == "Slightly Unique")/21819)*100),
                      ((sum(results_mar$Category == "Strongly Unique")/21819)*100))
)

#Cell_distribution

#Figure 1: Analysis geographical overview ----
##Panel A: Species richness ----
if(!"Phylacine_rast_Poj_Binary_SR.tif" %in% dir()){
  #Only living terrestrial and marine species included
  All_species <- c(Terrestrial_species, Marine_species) 
  Phylacine_rast_Poj_Binary_SR <- subset(Phylacine_rast_all_Poj_Binary, All_species)
  writeRaster(Phylacine_rast_Poj_Binary_SR,
              "Phylacine_rast_Poj_Binary_SR.tif",
              overwrite=T)
  
} else {
  Phylacine_rast_Poj_Binary_SR <-rast("Phylacine_rast_Poj_Binary_SR.tif")
}

#Sum all layers Phylacine_rast_Poj_Binary_SR to create SR raster
Species_richness <- sum(Phylacine_rast_Poj_Binary_SR)

#Create DF from raster
SR_df <- as.data.frame(Species_richness, xy = TRUE)
colnames(SR_df) <- c("lon", "lat", "SR")

#Align world map with raster
world_SR <- st_transform(world, crs(Species_richness))

#Plot SR DF
F1_A <-ggplot() +
  geom_raster(
    data = SR_df,
    aes(x = lon, y = lat, fill = SR)
  ) +
  geom_sf(
    data = world_SR,
    fill = NA,
    color = "grey",
    alpha = 1,
    linewidth = 0.09
  ) +
  coord_sf(expand = FALSE) +
  scale_fill_viridis_c(
    name = "Species richness",
    option = "D",
    limits = c(0, NA)
  ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold")) +
  labs(
    title = "(a)",
    x = NULL,
    y = NULL
  )

##Panel B: EN/CR species richness ----
if(!"Phylacine_rast_Poj_Binary_EN_SR.tif" %in% dir()){
  #Only EN/CR terrestrial and marine species included
  above_EN <- c(above_EN_Terr, above_EN_Mar)
  Phylacine_rast_Poj_Binary_EN_SR <- subset(Phylacine_rast_all_Poj_Binary, above_EN)
  writeRaster(Phylacine_rast_Poj_Binary_EN_SR,
              "Phylacine_rast_Poj_Binary_EN_SR.tif",
              overwrite=T)
  
} else {
  Phylacine_rast_Poj_Binary_EN_SR <-rast("Phylacine_rast_Poj_Binary_EN_SR.tif")
}

#Sum all layers Phylacine_rast_Poj_Binary_EN_SR to create SR raster
EN_species_richness <- sum(Phylacine_rast_Poj_Binary_EN_SR)

#Create DF from raster
EN_SR_df <- as.data.frame(EN_species_richness, xy = TRUE)
colnames(EN_SR_df) <- c("lon", "lat", "SR")

#Align world map with raster
world_EN_SR <- st_transform(world, crs(EN_species_richness))

#Plot SR DF
F1_B <-ggplot() +
  geom_raster(
    data = EN_SR_df,
    aes(x = lon, y = lat, fill = SR)
  ) +
  geom_sf(
    data = world_EN_SR,
    fill = NA,
    color = "grey",
    alpha = 1,
    linewidth = 0.09
  ) +
  coord_sf(expand = FALSE) +
  scale_fill_viridis_c(
    name = "EN/CR \nspecies richness",
    option = "D",
    limits = c(0, NA)
  ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold")) +
  labs(
    title = "(b)",
    x = NULL,
    y = NULL
  )

##Panel C: Cells to use ----
# plot(To_Use_cells_Terr, main = "Cells to keep 'Terrestrial'")
# plot(To_Use_cells_Mar, main = "Cells to keep 'Marine'")

#Convert logic rasters (True/false -> 1/NA)
To_Use_cells_Terr_plot <- ifel(To_Use_cells_Terr, 1, NA)
To_Use_cells_Mar_plot <- ifel(To_Use_cells_Mar, 1, NA)

To_Use_cells_All <- cover(To_Use_cells_Terr_plot, To_Use_cells_Mar_plot)

#Create DF from raster 
To_Use_Terr <- as.data.frame(To_Use_cells_Terr_plot, xy = TRUE, na.rm = TRUE)
colnames(To_Use_Terr) <- c("lon", "lat", "Use")
To_Use_Terr$type <- "Terrestrial included"

To_Use_Mar <- as.data.frame(To_Use_cells_Mar_plot, xy = TRUE, na.rm = TRUE)
colnames(To_Use_Mar) <- c("lon", "lat", "Use")
To_Use_Mar$type  <- "Marine included"

To_Use_All <- rbind(To_Use_Terr, To_Use_Mar)

#Align world map
world_TU <- st_transform(world, crs(To_Use_cells_All))

#Plot To_Use DF
F1_C <- ggplot() +
  geom_raster(
    data = To_Use_All,
    aes(x = lon, y = lat, fill = type)
  ) +
  geom_sf(
    data = world_TU,
    fill = NA,
    color = "grey50",
    alpha = 0.2,
    linewidth = 0.15
  ) +
  coord_sf(expand = FALSE) +
  scale_fill_manual(
    name = "Grid cells included",
    values = c(
      "Terrestrial included" = "#6CCE59",  # viridis green
      "Marine included"      = "#482878"   # viridis blue
    )
  ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold")) +
  labs(
    title = "(c)",
    x = NULL,
    y = NULL
  )

###Combine Panels ----
Figure1 <- F1_A / F1_B /F1_C
#plot(Figure1)


#Figure 2-3(A): Heat map ----
##Figure 2A: Terrestrial ----
F2_A <- ggplot(results_terr, aes(x = SES.Mean, y = SES.SD)) +
  geom_hex(bins = 30) +  # Hexagonal binning for density (Each hexagon represents amount of cells with similar values)
  scale_fill_viridis_c(name = "Cell Count",
                       limits = c(0, NA)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red2") + #Red line at 0
  labs(
    title = "(a)",
    x = "SES mean",
    y = "SES sd",) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))

##Figure 3A: Marine ----
F3_A <-ggplot(results_mar, aes(x = SES.Mean, y = SES.SD)) +
  geom_hex(bins = 30) +  # Hexagonal binning for density (Each hexagon represents amount of cells with similar values)
  scale_fill_viridis_c(name = "Cell Count",
                       limits = c(0, NA)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red2") + #Red line at 0
  labs(
    title = "(a)",
    x = "SES mean",
    y = "SES sd",) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))


#Figure 2-3 (B): Choropleth maps ----
##Empty global raster
tempwrlMap <-rast(ncols=180*2, nrows=180) #creates a blank raster grid covering the globe.
#(ncols=180*2, nrows=180) → grid resolution (here: 360 x 180 = 1° resolution)
tempwrlMap[] <- NA #fills all the cells with the value 'NA'

##Figure 2B: Terrestrial ----
#Extract cell numbers from Cell.ID in result dataframe
Cell_numbers_terr <- as.numeric(gsub("Cell.ID.","", results_terr$Cell.ID))

#Prepare a vector of NA values (for all raster cells)
vals_terr <- rep(NA, terra::ncell(tempwrlMap))

#Assign SES.mean values to only designated cells
vals_terr[Cell_numbers_terr] <- results_terr$SES.Mean

#convert dataframe for ggplot
SES_raster_terr <- setValues(tempwrlMap, vals_terr)
names(SES_raster_terr) <- "SES.mean"

###Map with coastlines
#range(results_terr$SES.Mean) #check SES range (-2.071278 to 5.485396)
F2_B <- ggplot() + 
  geom_spatraster(data = SES_raster_terr) +
  #Add country boundaries
  geom_sf(data = world, fill = NA, 
          color = "grey80",
          alpha = 0.2, 
          linewidth = 0.15) +
  #Create scale legend
  scale_fill_gradient2( 
    low = "#482878",
    mid = "grey95",
    high = "#6CCE59",
    midpoint = 0,
    na.value = "white",
    name = "SES.mean",
    limits = c(-3, 5.5)
  ) +
  #Overlay world map
  coord_sf(expand = FALSE) +
  labs(
    title = "(b)"
  ) +
  theme_minimal() + 
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right",
    panel.grid = element_line(color = "grey90", linewidth = 0.15),
    panel.background = element_rect(fill = NA, color = NA),
    panel.ontop = TRUE
  )

##Figure 3B: Marine ----
#Extract cell numbers from Cell.ID in result dataframe
Cell_numbers_mar <- as.numeric(gsub("Cell.ID.","", results_mar$Cell.ID))

#Prepare a vector of NA values (for all raster cells)
vals_mar <- rep(NA, terra::ncell(tempwrlMap))

#Assign SES.mean values to only designated cells
vals_mar[Cell_numbers_mar] <- results_mar$SES.Mean

#convert dataframe to spat raster for ggplot
SES_raster_mar <- setValues(tempwrlMap, vals_mar)
names(SES_raster_mar) <- "SES.mean"

###Map with coastlines
#range(results_mar$SES.Mean) #check SES range (-1.926227 to 9.620645)
F3_B <- ggplot() + 
  geom_spatraster(data = SES_raster_mar) +
  #Add country boundaries
  geom_sf(data = world, fill = NA, color = "grey80",alpha = 0.2, linewidth = 0.15) +
  scale_fill_gradient2(
    low = "#482878",
    mid = "grey95",
    high = "#6CCE59",
    midpoint = 0,
    na.value = "white",
    name = "SES.mean",
    limits = c(-4, 10)
  ) +
  coord_sf(expand = FALSE) +
  labs(
    title = "(b)"
  ) +
  theme_minimal() + 
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right",
    panel.grid = element_line(color = "grey90", linewidth = 0.15),
    panel.background = element_rect(fill = NA, color = NA),
    panel.ontop = TRUE
  )

#Combine Panels figure 2 & 3 ----
#Over each other
# Figure2 <- F2_A / F2_B 
# Figure3 <- F3_A / F3_B
#Side by side
Figure2 <- F2_A | F2_B
Figure3 <- F3_A | F3_B
#plot(Figure2)
#plot(Figure3)


#Figure 4: Extreme cells ----
##Panel A + B ----
#Find Extreme cells
  #Terrestrial DF
  #Create dataframe containing extreme SES cells
extreme_ses_terr <- results_terr[abs(results_terr$SES.Mean) > 5, ] #abs returns absolute value (removes -, catches both extreme positives and negatives)
#Print overlook
# if(nrow(extreme_ses_terr) > 0) {
#   print(extreme_ses_terr[, c("Cell.ID", "Species.Richness", "EN.CR.Richness", "SES.Mean", "SES.SD")])
# }

  #Marine DF
  #Create dataframe containing extreme SES cells
extreme_ses_mar <- results_mar[abs(results_mar$SES.Mean) > 5, ] #abs returns absolute value (removes -, catches both extreme positives and negatives)
#Print overlook
# if(nrow(extreme_ses_mar) > 0) {
#   print(extreme_ses_mar[, c("Cell.ID", "Species.Richness", "EN.CR.Richness", "SES.Mean", "SES.SD")])
# }

  #Create raster
  #Terrestrial raster
#Extract cell numbers from Cell.ID
Cell_numbers_EX_terr <- as.numeric(gsub("Cell.ID.","", extreme_ses_terr$Cell.ID))

#Prepare a vector of NA values (for all raster cells)
vals_EX_terr <- rep(NA, terra::ncell(tempwrlMap))

#Assign value to only designated cells
#vals_EX_terr[Cell_numbers_EX_terr] <- extreme_ses_terr$SES.Mean
vals_EX_terr[Cell_numbers_EX_terr] <- 1  # Binary: 1 = extreme cell

#Convert dataframe to spat raster for ggplot
SES_raster_EX_terr <- setValues(tempwrlMap, vals_EX_terr)
#names(SES_raster_EX_terr) <- "SES.mean"
names(SES_raster_EX_terr) <- "Extreme"

#SES_raster_EX_terr_Binary <- (SES_raster_EX_terr>0)*1
#plot(sum(SES_raster_EX_terr))


  #Marine raster
#Extract cell numbers from Cell.ID
Cell_numbers_EX_mar <- as.numeric(gsub("Cell.ID.", "", extreme_ses_mar$Cell.ID))

#Prepare a vector of NA values (for all raster cells)
vals_EX_mar <- rep(NA, terra::ncell(tempwrlMap))

#Assign value to only designated cells
vals_EX_mar[Cell_numbers_EX_mar] <- 1  # Binary: 2 = extreme cell (different value for distinction)

#Convert dataframe to spat raster for ggplot
SES_raster_EX_mar <- setValues(tempwrlMap, vals_EX_mar)
names(SES_raster_EX_mar) <- "Extreme"

  #Combined raster
# Combine rasters (marine will overlay terrestrial where both exist)
combined_raster <- SES_raster_EX_terr
#combined_raster <- SES_raster_EX_mar
combined_raster[!is.na(SES_raster_EX_mar)] <- SES_raster_EX_mar[!is.na(SES_raster_EX_mar)]

combined_raster <- as.factor(combined_raster)

#Plot
###A - Terrestrial extreme cells (Southeast Asia focus) ----
SES_raster_EX_terr_f <- as.factor(SES_raster_EX_terr)

F4A <- ggplot() +
  geom_spatraster(data = SES_raster_EX_terr_f) +
  #Green color for terrestrial extreme cells
  scale_fill_manual(
    values = c("1" = "#6CCE59"), 
    labels = c("Terrestrial extreme cells"),
    na.translate = FALSE, #Don't show NA values in legend
    name = "" #Remove legend title
  ) +
  #Add country boundaries
  geom_sf(data = world,
          fill = NA, color = "grey70", linewidth = 0.2) +
  #Zoom to Southeast Asia region
  coord_sf(xlim = c(80, 135),   #Longitude range
           ylim = c(-15, 35),   #Latitude range
           expand = FALSE) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold")
  ) +
  labs(title = "(a)")

###B - Marine extreme cells (Galápagos/Eastern Pacific focus)----
SES_raster_EX_mar_f <- as.factor(SES_raster_EX_mar)

F4B <- ggplot() +
  geom_spatraster(data = SES_raster_EX_mar_f) +
  #Purple color for marine extreme cells
  scale_fill_manual(
    values = c("1" = "#482878"), 
    labels = c("Marine extreme cells"),
    na.translate = FALSE, #Don't show NA values in legend
    name = "" #Remove legend title 
    ) +
  #Add country boundaries
  geom_sf(data = world,
          fill = NA, color = "grey70", linewidth = 0.2) +
  #Zoom to Galápagos/Eastern Pacific region
  coord_sf(xlim = c(-115, -65),  #Longitude range
           ylim = c(-25, 25),    #Latitude range
           expand = FALSE) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold")
  ) +
  labs(title = "(b)")

#Combine Facets
#Side by side
F4_maps <- F4A | F4B 
#plot(F4_maps)


##Panel C: Terrestrial species composition ----

###Species in-common ----
#This creates a phylogeny of species that occur in EVERY extreme cell (intersection)

#Get species present in ALL extreme terrestrial cells (intersection)
all_species_terr <- Reduce(intersect, lapply(extreme_ses_terr$Cell.ID, function(cell_id) {
  cell_idx <- which(names(Spp_per_Use_cells_Terr) == cell_id)
  Spp_per_Use_cells_Terr[[cell_idx]]$Spp
}))
#reduce(intersect, list) finds species present in ALL lists (intersection)
#unique() keeps only species common to all cells

#Get ALL EN/CR species from extreme cells
all_encr_terr <- unique(unlist(lapply(extreme_ses_terr$Cell.ID, function(cell_id) {
  cell_idx <- which(names(Spp_per_Use_cells_Terr) == cell_id)
  cell_data <- Spp_per_Use_cells_Terr[[cell_idx]]
  cell_data$Spp[cell_data$IUCN %in% c("EN", "CR")]
})))

#Filter EN/CR to only those present in all cells
all_encr_terr <- all_encr_terr[all_encr_terr %in% all_species_terr]

#Prune multiphylo (1000 trees) to species present in ALL extreme cells
pruned_trees_terr <- list()
for (d in 1:length(Complete_trees)) {
  if(d %% 200 == 0) cat("  Pruning tree", d, "/", length(Complete_trees), "\n")
  
  current_tree <- Complete_trees[[d]]
  species_to_keep <- all_species_terr[all_species_terr %in% current_tree$tip.label]
  
  #Only prune if we have more than 1 species to keep
  if (length(species_to_keep) > 1) {
    pruned_tree <- drop.tip(current_tree, 
                            current_tree$tip.label[!current_tree$tip.label %in% species_to_keep])
    pruned_trees_terr[[d]] <- pruned_tree
  }
}

#Convert list to multiPhylo object for consensus tree calculation
class(pruned_trees_terr) <- "multiPhylo"

#Create 50% majority-rule consensus tree
  #p = 0.5: branches appearing in >50% of trees are retained
consensus_terr <- consensus(pruned_trees_terr, p = 0.5, rooted = TRUE)

#Save consensus tree
write.nexus(consensus_terr, file = "Consensus_Tree_Terrestrial_Extreme_AllCells.nex")

#Plot consensus tree
png("Terrestrial_Extreme_Consensus_Tree_AllCells.png", width = 1800, height = 1800, res = 150) #Open graphics device

#Color EN/CR species in red, others in black
tip_colors <- rep("black", length(consensus_terr$tip.label))
tip_colors[consensus_terr$tip.label %in% all_encr_terr] <- "red"

#Plot as fan tree
plot(consensus_terr, 
     type = "fan",
     tip.color = tip_colors,  #Species color coding
     cex = 1.0)

#Add legend explaining color coding
legend("topleft", 
       legend = c(paste("EN/CR species (n=", length(all_encr_terr), ")"),
                  paste("Other species (n=", length(all_species_terr) - length(all_encr_terr), ")")),
       col = c("red", "black"),
       lty = 1, lwd = 2,
       cex = 1.0)

dev.off() #Close the graphics device

#Read back the saved tree
consensus_terr_ALL <-read.nexus("Consensus_Tree_Terrestrial_Extreme_AllCells.nex")


##Panel D: Marine species composition ----

###Species in-common ----
#if species is present in cell 1 AND cell 2 AND cell 3 AND ...

#Get species present in ALL extreme marine cells (intersection)
all_species_mar <- Reduce(intersect, lapply(extreme_ses_mar$Cell.ID, function(cell_id) {
  cell_idx <- which(names(Spp_per_Use_cells_Mar) == cell_id)
  Spp_per_Use_cells_Mar[[cell_idx]]$Spp
}))
#reduce(intersect, list) finds species present in ALL lists (intersection)
#unique() keeps only species common to all cells

#Get ALL EN/CR species from extreme cells
all_encr_mar <- unique(unlist(lapply(extreme_ses_mar$Cell.ID, function(cell_id) {
  cell_idx <- which(names(Spp_per_Use_cells_Mar) == cell_id)
  cell_data <- Spp_per_Use_cells_Mar[[cell_idx]]
  cell_data$Spp[cell_data$IUCN %in% c("EN", "CR")]
})))

#Filter EN/CR to only those present in all cells
all_encr_mar <- all_encr_mar[all_encr_mar %in% all_species_mar]

#Prune multiphylo (1000 trees) to species present in ALL extreme cells
pruned_trees_mar <- list()
for (d in 1:length(Complete_trees)) {
  if(d %% 200 == 0) cat("  Pruning tree", d, "/", length(Complete_trees), "\n")
  
  current_tree <- Complete_trees[[d]]
  species_to_keep <- all_species_mar[all_species_mar %in% current_tree$tip.label]
  
  #Only prune if we have more than 1 species to keep
  if (length(species_to_keep) > 1) {
    pruned_tree <- drop.tip(current_tree, 
                            current_tree$tip.label[!current_tree$tip.label %in% species_to_keep])
    pruned_trees_mar[[d]] <- pruned_tree
  }
}

#Convert list to multiPhylo object for consensus tree calculation
class(pruned_trees_mar) <- "multiPhylo"

#Create 50% majority-rule consensus tree
  #p = 0.5: branches appearing in >50% of trees are retained
consensus_mar <- consensus(pruned_trees_mar, p = 0.5, rooted = TRUE)

#Save consensus tree
write.nexus(consensus_mar, file = "Consensus_Tree_Marine_Extreme_AllCells.nex")

#Plot consensus tree
png("Marine_Extreme_Consensus_Tree_AllCells.png", width = 1800, height = 1800, res = 150) #Open graphics device

#Color EN/CR species in red, others in black
tip_colors <- rep("black", length(consensus_mar$tip.label))
tip_colors[consensus_mar$tip.label %in% all_encr_mar] <- "red"

#Plot as fan tree
plot(consensus_mar, 
     type = "fan",
     tip.color = tip_colors,
     cex = 1.0)

#Add legend explaining color coding
legend("topleft", 
       legend = c(paste("EN/CR species (n=", length(all_encr_mar), ")"),
                  paste("Other species (n=", length(all_species_mar) - length(all_encr_mar), ")")),
       col = c("red", "black"),
       lty = 1, lwd = 2,
       cex = 1.0)

dev.off() #Close the graphics device

#Read back the saved tree
consensus_mar_ALL <-read.nexus("Consensus_Tree_Marine_Extreme_AllCells.nex")



#Plot Figures ----
#Figure 1 - Species distribution Maps
plot(Figure1)

#Figure 2 - Terrestrial
plot(Figure2)

#Figure 3 - Marine
plot(Figure3)

#Figure 4 - Extreme cells
plot(F4_maps)
