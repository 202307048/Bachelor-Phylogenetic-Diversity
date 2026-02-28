#Analysis of Terrestrial Mammals

#Terrestrial 'Raster'----
#Creates a binary and projected raster stack containing the Terrestrial species
if(!"Phylacine_rast_Terrestrial_Poj_Binary.tif" %in% dir()){
  Phylacine_rast_Terrestrial_Poj_Binary <- subset(Phylacine_rast_all_Poj_Binary, Terrestrial_species)
  plot(sum(Phylacine_rast_Terrestrial_Poj_Binary), main = "Terrestrial")
  writeRaster(Phylacine_rast_Terrestrial_Poj_Binary,
              "Phylacine_rast_Terrestrial_Poj_Binary.tif",
              overwrite=T)
  
} else {
  Phylacine_rast_Terrestrial_Poj_Binary <-rast("Phylacine_rast_Terrestrial_Poj_Binary.tif")
}

##Above EN Terrestrial 'Raster'----
#Creates a binary and projected raster stack containing the Terrestrial species with a IUCN rating above EN
if(!"Phylacine_rast_Terrestrial_EN_Poj_Binary.tif" %in% dir()){
  Phylacine_rast_Terrestrial_EN_Poj_Binary <- subset(Phylacine_rast_all_Poj_Binary, above_EN_Terr)
  plot(sum(Phylacine_rast_Terrestrial_EN_Poj_Binary), main = "Above EN Terrestrial")  
  writeRaster(Phylacine_rast_Terrestrial_EN_Poj_Binary,
            "Phylacine_rast_Terrestrial_EN_Poj_Binary.tif",
            overwrite =T)
} else {
  Phylacine_rast_Terrestrial_EN_Poj_Binary <-rast("Phylacine_rast_Terrestrial_EN_Poj_Binary.tif")
}

#Find usable grid cells - global raster map ----
#plot(sum(Phylacine_rast_Terrestrial_Poj_Binary), main ="Terrestrial"), adds up the raster → produces a species richness map (number of terrestrial species present per cell), plots the result
SR_Terr <- sum(Phylacine_rast_Terrestrial_Poj_Binary) #Terrestrial SR sum for each grid

##Subset grids with ≤ 10 species richness ----
low_SR_Terr <- SR_Terr <= 10 #identify cells with 10 or fewer SR
#plot(low_SR_Terr, main = "Cells with ≤10 Terrestrial species") #visualize distribution of low SR grids

##Subset grids with ≤ 2 endangered species ----
EN_richness_Terr <- sum(Phylacine_rast_Terrestrial_EN_Poj_Binary) #adds up endangered species layers to get endangered species richness per cell

low_EN_Terr <- EN_richness_Terr <= 2 #identify cells with 2 or fewer endangered species
#plot(low_EN_Terr, main = "Cells with ≤2 endangered Terrestrial species") #visualize distribution of low EN grids

###Merge low-richness and low-endangerment cells, remove from Global raster map ----
To_Use_cells_Terr <- (low_SR_Terr != 1) & (low_EN_Terr != 1) #creates list excluding low SR/EN grids
# Create a logical mask that identifies cells to keep:
# TRUE = cell passes thresholds (SR >10 AND EN >2)
# FALSE = cell excluded
#plot(To_Use_cells_Terr, main = "Cells to keep 'Terrestrial'") #visualize distribution of 'valid' cells that remain 

#names(values(Phylacine_rast_Terrestrial_EN_Poj_Binary)[i,])[values(Phylacine_rast_Terrestrial_EN_Poj_Binary)[i,]!=0]
  #Uses the TRUE/FALSE vector to subset the species names
  #Returns ONLY the names where value ≠ 0 (i.e., species present)

##Extract ID's of usable/'valid' cells ----
Cells_to_use_Terr <- which(values(To_Use_cells_Terr)==T) #converts logical raster into a list of cell indices (numeric IDs) 

##Extract species composition for each grid cell ----
Phylacine_DataFrame_Terrestrial_Poj_Binary <- values(Phylacine_rast_Terrestrial_Poj_Binary) #converts raster stack into dataframe

##Spp_per_Use_cells_Terr ----
#Loop through each "usable/valid" grid cell
#Extract which species are present (value = 1)
#Match each species to its IUCN status from Trait_data_filtered
Spp_per_Use_cells_Terr <- lapply(Cells_to_use_Terr,
                            function(i){#i<-Cells_to_use_Terr[1]
                              Spp_List <- Phylacine_DataFrame_Terrestrial_Poj_Binary[i,]
                              #Extract species presence/absence for this cell
                              Out <- data.frame(Spp=names(Spp_List[Spp_List!=0]))
                              #Keep only species that are present (value ≠ 0)
                              Out$IUCN <- Trait_data_filtered$IUCN.Status.1.2[match(Out$Spp,Trait_data_filtered$Binomial.1.2)]
                              #Match IUCN status using Trait_data_filtered and add as new column
                              return(Out)
                              #Return a dataframe: species name + IUCN category
                            })
names(Spp_per_Use_cells_Terr) <- paste0("Cell.ID.",Cells_to_use_Terr) #Name each list element by its corresponding cell ID
Spp_per_Use_cells_Terr[1:10] #Inspect the first 10 elements

#Output - Data frame (Spp_per_Use_cells_Terr)
  #List of data frames (cell ID's) 
    #For each cell ID a data frame: species name + matching IUCN category
      #$Spp - vector of species present in the grid
      #$IUCN - IUCN assessment of each species ('LC', 'VU', 'EN', 'CR')


#ANALYSIS Terrestrial ----

##Goal list ----
  # Prune multiphylo to species present in grid cell
  # Estimate PD - full 
  # Estimate PD - minus EN species
  # Estimate PD change (full - EN)
  # Random communities (1000 communities)
    # Remove amount of species, corresponding to amount of EN species present, at rand 
    # Estimate PD -  minus Rnd species removal (Rnd_EN)
  # Estimate PD rnd-change (full - Rnd-EN removal)
    # done 1000 times (Once for each community)
  # Estimate Rnd Mean of PD rnd-change
  # Estimate Rnd SD of PD rnd-change
  # Estimate SES
  # Save as table (.csv)
    # cellID, PD_full, PD_noEN, MEAN_Rnd_EN_Rem, SD_Rnd_EN_Rem, SES_mean, SES_sd

# Load required packages
#library(terra)         #in Basis script
#library(picante)       #in Basis script
#library(ape)           #in Basis script
#library(parallel)      #in Basis script
#library(foreach)       #in Basis script
#library(doParallel)    #in Basis script

#Spp_per_Use_cells_Terr = list of data frames, each element corresponds to one of the filtered grid cells (≤ 10 species + ≤ 2 EN species)
#Complete_trees = multiphylo object (list of 1000 full phylogenetic trees)


## Optimized function for phylogenetic diversity analysis ----
  #Creates analysis function
calculate_ses_pd <- function(cell_data, tree_list, n_rand = 1000, 
                             output_file = "SES_Summ_Final.csv",
                             n_cores = NULL, chunk_size = 100,
                             prune_trees = TRUE) {
  
  #Setup parallel backend 
  if (is.null(n_cores)) {
    #If number of cores aren't specified it detects number of available cores minus 1 (keeping computer responsive)
    n_cores <- max(1, detectCores() - 1)
  }
  #Create cluster of worker processes (one pr. core)
  cl <- makeCluster(n_cores)
  #Register cluster so foreach knows to use it
  registerDoParallel(cl)
  
  #Export required packages to workers
  clusterEvalQ(cl, {
    library(picante) #Needed for pd() function
    library(ape)     #Needed for drop.tip() function (pruning)
  })
  
  #Print analysis parameters to console (Not needed - excess)
  cat(sprintf("Processing %d cells using %d cores with %d trees and %d random communities\n", 
              length(cell_data), n_cores, length(tree_list), n_rand))
  cat("================================================================================\n") #create section divide to improve 'readability' when printed in console
  
  #Initialize output file
    #If file exists from previous run delete it and start fresh
  if (file.exists(output_file)) file.remove(output_file)
  
  #Process cells in chunks to manage memory
  n_cells <- length(cell_data) #Total number of grid cells (2514)
  n_chunks <- ceiling(n_cells / chunk_size) #Calculate amount of chunks needed
  overall_start <- Sys.time() #Track total processing time
  
  #Loop through chunks (e.g. cells 1-50 or 51-100...)
  for (chunk_idx in 1:n_chunks) {
    #Calculate which cells are in this chunk
    chunk_start <- (chunk_idx - 1) * chunk_size + 1 #First cell in chunk
    chunk_end <- min(chunk_idx * chunk_size, n_cells) #Last cell in chunk
    chunk_indices <- chunk_start:chunk_end #Creates vector of cell indices to process
    
    chunk_time_start <- Sys.time() #Track cell processing time
    
    #Process chunk in parallel
      #Foreach distributes cells (chunks) across worker cores
    chunk_results <- foreach(i = chunk_indices, #Each worker gets different values of i (different chunks of cell indices to process)
                             .combine = rbind, #Combines results from workers by rowbinding
                             .packages = c("picante", "ape"), #Packages needed by workers
                             .errorhandling = "pass") %dopar% { #"pass" = continue if one cell fails
                                                                #"%dopar%" = run the following in parallel: 
                               
                               #Wrapped in trycatch to handle error 'gracefully'
                               tryCatch({
                                 cell_start <- Sys.time() #Track time for this cell
                                 Comm <- cell_data[[i]] #Extract species assemblage for this cell
                                 
                                 #Calculate species richness metrics
                                 total_richness <- nrow(Comm) #Species richness in the cell
                                 en_cr_richness <- sum(Comm$IUCN %in% c("EN", "CR")) #Number of Endangered species in this cell
                                 
                                 #Setup community matrix
                                  #Create presence/absence columns for different scenarios (original species assemblage, species assemblage without endangered)
                                 Comm$Org <- 1 #'Original' - all species present (value = 1 - present)
                                 Comm$No.EN <- 1 #No endangered (EN & CR) community (value = 1 - present)
                                 Comm$No.EN[Comm$IUCN %in% c("EN", "CR")] <- 0 #Remove EN/CR species (sets their value to 0 - absent) from the 'No endangered' community
                                 
                                 #Pre-allocate random community matrix
                                 n_species <- nrow(Comm) #Number of species in cell
                                 n_to_remove <- en_cr_richness #Number of species to remove in each random community
                                 #Create Matrix: rows = species, columns = random communities (1-1000), all values = 1
                                 rand_matrix <- matrix(1, nrow = n_species, ncol = n_rand) #n_rand - number of random communities wanted
                                 
                                 #Vectorized random community generation
                                  #For each random community, remove the same number of species, as number of EN/CR present, but randomly
                                 if (n_to_remove > 0) {
                                   for (j in 1:n_rand) { #Loop through each random community
                                     #Randomly select which species to remove
                                     to_remove <- sample.int(n_species, n_to_remove)
                                     #Set choosen species value to 0 (absent) in this random community
                                     rand_matrix[to_remove, j] <- 0
                                   }
                                 }
                                 
                                 #Combine original and random communities
                                  #Columns: Org (all present), No.EN (EN/CR removed), Rand.1, Rand.2, ..., Rand.1000
                                 comm_matrix <- cbind(Org = Comm$Org, No.EN = Comm$No.EN, rand_matrix)
                                 
                                 #Ensure species names are unique and valid (Extra safety - double check)
                                 spp_names <- as.character(Comm$Spp) #Get species names as characters 
                                 if(any(duplicated(spp_names))) { #Check for duplicates - (0 found)
                                   #If duplicates exist, make unique by adding .1, .2, etc.
                                   spp_names <- make.unique(spp_names) #
                                 }
                                 rownames(comm_matrix) <- spp_names #set species as rownames in matrix
                                 #Transpose matrix for pd() function (needs communities as rows and species as columns)
                                 comm_matrix_t <- t(comm_matrix) 
                                 
                                 #Ensure community names are unique for pd() function
                                 rownames(comm_matrix_t) <- make.unique(c("Org", "No.EN", paste0("Rand.", 1:n_rand)))
                                 
                                 #Calculate PD for all trees
                                  #Pre-allocate vectors to store results from each tree
                                 n_trees <- length(tree_list) #Number of trees in multiphylo
                                  #Create empty numeric vectors
                                 ses_values <- numeric(n_trees) 
                                 pd_org_all <- numeric(n_trees)
                                 pd_no_en_all <- numeric(n_trees)
                                 mean_pd_rnd_all <- numeric(n_trees)
                                 sd_pd_rnd_all <- numeric(n_trees)
                                 
                                 #Loop through each phylogenetic tree (d is the tree index)
                                  #Prune each tree to only species in this community
                                 for (d in 1:n_trees) { 
                                   current_tree <- tree_list[[d]]
                                   if (prune_trees) {
                                     #Find which species in this cell exist in the current tree
                                     species_to_keep <- spp_names[spp_names %in% current_tree$tip.label]
                                     if (length(species_to_keep) > 1) { #(Extra safety - double check species richness in cell)
                                       #Remove all species not in this cell's community
                                       current_tree <- drop.tip(current_tree, 
                                                                current_tree$tip.label[!current_tree$tip.label %in% species_to_keep])
                                     }
                                   }
                                   
                                   #Calculate PD
                                    #calculate PD for each community in matrix
                                      #include.root = FALSE - measures from MRCA of present species, not tree root
                                   pd_vals <- pd(comm_matrix_t, current_tree, include.root = FALSE)
                                   
                                   #Extract values
                                    #pd_vals dataframe contain rows for each community and columns $PD and $SR
                                   pd_org <- pd_vals$PD[1] #PD of row 1 (Original community)
                                   pd_no_en <- pd_vals$PD[2] #PD of row 2 (No EN/CR community)
                                   pd_rand <- pd_vals$PD[3:nrow(pd_vals)] #PD of rows 3 onward (all random communities)pd_rand - vector of PD results)
                                   
                                   #Calculate losses
                                    #How much PD is lost when species are removed
                                   observed_loss <- pd_org - pd_no_en  #loss from removing EN/CR
                                   random_losses <- pd_org - pd_rand   #losses from removing random species (vector)
                                   
                                   #Store PD metrics
                                   pd_org_all[d] <- pd_org #Store PD of original community
                                   pd_no_en_all[d] <- pd_no_en #Store PD without EN/CR
                                   mean_pd_rnd_all[d] <- mean(pd_rand) #Store average PD of random communities
                                   sd_pd_rnd_all[d] <- sd(pd_rand) #Store variability of random community PD (standard deviation)
                                   
                                   #Calculate SES (how unusual is 'observed loss' vs random losses)
                                   mean_random_loss <- mean(random_losses) #Average PD loss from removing random species
                                   sd_random_loss <- sd(random_losses) #Standard deviation of random losses
                                   ses_values[d] <- (observed_loss - mean_random_loss) / sd_random_loss #SES
                                 } #End of tree loop -  results from all trees
                                 
                                 #Summarize results
                                  #Average the results across all trees (accounts for phylogenetic uncertainty)
                                 result <- data.frame(
                                   Cell.ID = names(cell_data)[i],       #Cell identifier (ID)
                                   Species.Richness = total_richness,   #Total number of species in cell
                                   EN.CR.Richness = en_cr_richness,     #Number of EN/CR species in cell
                                   SES.Mean = mean(ses_values),         #Average SES across all trees (Main result)
                                   SES.SD = sd(ses_values),             #How much SES varies across trees
                                   PD.Org.Mean = mean(pd_org_all),      #Average PD of original community 
                                   PD.Org.SD = sd(pd_org_all),          #Variability in original PD across trees
                                   PD.NoEN.Mean = mean(pd_no_en_all),   #Average PD without EN/CR
                                   PD.NoEN.SD = sd(pd_no_en_all),       #Variability in PD without EN/CR
                                   PD.Rnd.Mean = mean(mean_pd_rnd_all), #Average of random community PD means
                                   PD.Rnd.SD = mean(sd_pd_rnd_all),     #Average of random community PD SDs
                                   stringsAsFactors = FALSE #Don't convert strings to factors
                                 )
                                 
                                 #Add processing time as attribute (not in output)
                                 attr(result, "proc_time") <- as.numeric(Sys.time() - cell_start, units = "secs")
                                 
                                 #Memory clean up
                                  #remove large objects to free memory before processing next cell
                                 rm(comm_matrix, comm_matrix_t, ses_values, pd_org_all, 
                                    pd_no_en_all, mean_pd_rnd_all, sd_pd_rnd_all)
                                  #Force garbage collection to free memory immediately
                                 gc(verbose = FALSE) #verbose = FALSE - don't print memory stats to console
                                 
                                 #return results data frame for this cell
                                 result
                                 
                               }, #End of tryCatch block
                               error = function(e) { #(Extra safety)
                                 #If anything goes wrong in the try block above, this code runs instead
                                 #Returns error information for debugging
                                 data.frame(
                                   Cell.ID = names(cell_data)[i],
                                   Species.Richness = NA,
                                   EN.CR.Richness = NA,
                                   SES.Mean = NA,
                                   SES.SD = NA,
                                   PD.Org.Mean = NA,
                                   PD.Org.SD = NA,
                                   PD.NoEN.Mean = NA,
                                   PD.NoEN.SD = NA,
                                   PD.Rnd.Mean = NA,
                                   PD.Rnd.SD = NA,
                                   #stores error message 
                                   Error = paste("ERROR in cell", i, ":", as.character(e)),
                                   stringsAsFactors = FALSE
                                 )
                               }) #End of error function and tryCatch
                             } #End of foreach loop - all cells in chunk processed
    
    #Write chunk results to file
    write.table(chunk_results, output_file, #Data frame of results from this chunk
                sep = ",", append = (chunk_idx > 1),
                  #if chunk 2+, append to existing file
                  #if chunk 1, create new file (overwrite)
                col.names = (chunk_idx == 1), row.names = FALSE) #Only write column headers for first chunk, don't write row numbers
    
    #Calculate and display progress
    chunk_time <- as.numeric(Sys.time() - chunk_time_start, units = "secs") #Time (seconds) it took to analyse chunk
    cells_processed <- chunk_end #Total number of cells completed so far
    cells_remaining <- n_cells - cells_processed #Number of cells left to process
    avg_time_per_cell <- chunk_time / nrow(chunk_results) #Average seconds per cell in this chunk
    estimated_remaining <- (cells_remaining / chunk_size) * chunk_time #Rough estimate of remaining time
    
    #Print progress bar to console
    cat(sprintf("Chunk %d/%d complete | Cells: %d-%d/%d (%.1f%%) | Time: %.1fs (%.2fs/cell) | Est. remaining: %.1f min\n", #%.1f - "number with 1 decimal place, %d - integer
                chunk_idx, n_chunks, chunk_start, chunk_end, n_cells,
                (cells_processed/n_cells)*100, chunk_time, avg_time_per_cell,
                estimated_remaining/60)) #estimate time left in minutes
    
    #Force garbage collection between chunks
    gc(verbose = FALSE)
    
  } #End of chunk loop - all chunks processed
  
  #Cleanup parallel cluster
    #Closes all worker processes and frees up CPU cores
  stopCluster(cl)
  
  #Calculate total processing time
  total_time <- as.numeric(Sys.time() - overall_start, units = "mins")
  
  #Display completion summary
  cat("================================================================================\n") #create section divide to improve 'readability' when printed in console
  cat(sprintf("Analysis complete Total time: %.2f minutes (%.2f hours)\n", 
              total_time, total_time/60)) #total time in hours
  cat(sprintf("Results saved to: %s\n", output_file)) #%s - placeholder for string (filename)
  
  #Return summary
  final_results <- read.csv(output_file) #Creates data frame from the finished file
  return(final_results)
    #Returns the data frame to the user
    #Allows immediately work with results without manually loading the CSV
  
} #End of function

## Function usage - workflow ----
if(!"Results_terr_df.rds"%in%dir()){ #Don't run if results already exist - avoid re-running analysis
  results_terr <- calculate_ses_pd(
    cell_data = Spp_per_Use_cells_Terr,         #Input: List of terrestrial cells
    tree_list = Complete_trees,                 #Input: 1000 phylogenetic trees
    n_rand = 1000,                              #1000 random communities per cell
    output_file = "SES_Summ_Final_Terr.csv",    #filename - save results to this file
    n_cores = 7,                                #Number of cores used for analysis
    chunk_size = 50,                            #Number of cells in chunks
    prune_trees = TRUE                          #Prune trees to cell assemblage
  )
  saveRDS(results_terr, file = "Results_terr_df.rds") #saves data frame as rds file (faster to load)
} else {
  results_terr <- readRDS("Results_terr_df.rds") #Load rds file if file already exist
}

## Output file structure ----
#The output CSV contains:
  #Cell.ID: Grid cell identifier
  #Species.Richness: Total number of species in cell
  #EN.CR.Richness: Number of EN/CR species in cell
  #SES.Mean: Mean standardized effect size 
  #SES.SD: Variability of SES across trees 
  #PD.Org.Mean: Mean PD of original community
  #PD.Org.SD: Variability in original PD
  #PD.NoEN.Mean: Mean PD without EN/CR species
  #PD.NoEN.SD: Variability in PD without EN/CR
  #PD.Rnd.Mean: Mean random community PD
  #PD.Rnd.SD: Mean variability in random PD

#Run history
# FIRST RUN
# Full data (1000 trees and 1000 Random communities using 6 cores)
# Output name:  "SES_Summ_Final_Terr.csv"
# Runtime:       43,8 hours
# Dates:         26.11.25 - 02 am
#                27.11.25 - 09 pm

# Error detected:
#   In the filtering of cells in Spp_per_Use_cells_Terr ('contaminated with marine species')
# - Had used Phylacine_rast_all_Poj_Binary instead of Phylacine_rast_Terrestrial_Poj_Binary
# Now corrected in both Terrestrial and Marine script

#FINAL RUN
# Full data (1000 trees and 1000 Random communities using 7 cores)
# Output name:  "SES_Summ_Final_Terr.csv"
# Runtime:      54,1 hours
# Dates:        27.11.25 - 12 am
#               30.11.25 - 6 am
