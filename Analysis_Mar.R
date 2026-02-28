#Analysis of Marine Mammals

#Marine 'Raster'----
#Creates a binary and projected raster stack containing the Marine species
if(!"Phylacine_rast_Marine_Poj_Binary.tif" %in% dir()){
  Phylacine_rast_Marine_Poj_Binary <- subset(Phylacine_rast_all_Poj_Binary, Marine_species)
  plot(sum(Phylacine_rast_Marine_Poj_Binary), main = "Marine")
  writeRaster(Phylacine_rast_Marine_Poj_Binary,
              "Phylacine_rast_Marine_Poj_Binary.tif",
              overwrite=T)
  
} else {
  Phylacine_rast_Marine_Poj_Binary <-rast("Phylacine_rast_Marine_Poj_Binary.tif")
}
#plot(sum(Phylacine_rast_Marine_Poj_Binary), main = "Marine")

##Above EN Marine 'Raster'----
#Creates a binary and projected raster stack containing the Marine species with a IUCN rating above EN
if(!"Phylacine_rast_Marine_EN_Poj_Binary.tif" %in% dir()){
  Phylacine_rast_Marine_EN_Poj_Binary <- subset(Phylacine_rast_all_Poj_Binary, above_EN_Mar)
  plot(sum(Phylacine_rast_Marine_EN_Poj_Binary), main = "Above EN Marine")  
  writeRaster(Phylacine_rast_Marine_EN_Poj_Binary,
              "Phylacine_rast_Marine_EN_Poj_Binary.tif",
              overwrite =T)
} else {
  Phylacine_rast_Marine_EN_Poj_Binary <-rast("Phylacine_rast_Marine_EN_Poj_Binary.tif")
}
#plot(sum(Phylacine_rast_Marine_EN_Poj_Binary), main = "Above EN Marine")


#Find usable grid cells - global raster map ----
#plot(sum(Phylacine_rast_Marine_Poj_Binary), main ="Marine"), adds up the raster → produces a species richness map (number of Marine species present per cell), plots the result
SR_Mar <- sum(Phylacine_rast_Marine_Poj_Binary) #Marine SR sum for each grid

##Subset grids with ≤ 10 species richness ----
low_SR_Mar <- SR_Mar <= 10 #identify cells with 10 or fewer SR
#plot(low_SR_Mar, main = "Cells with ≤10 Marine species") #visualize distribution of low SR grids

##Subset grids with ≤ 2 endangered species ----
EN_richness_Mar <- sum(Phylacine_rast_Marine_EN_Poj_Binary) #adds up endangered species layers to get endangered species richness per cell

low_EN_Mar <- EN_richness_Mar <= 2 #identify cells with 2 or fewer endangered species
#plot(low_EN_Mar, main = "Cells with ≤2 endangered Marine species") #visualize distribution of low EN grids

###Merge low-richness and low-endangerment cells, remove from Global raster map ----
To_Use_cells_Mar <- (low_SR_Mar != 1) & (low_EN_Mar != 1) #creates list excluding low SR/EN grids
# Create a logical mask that identifies cells to keep:
# TRUE = cell passes thresholds (SR >10 AND EN >2)
# FALSE = cell excluded
plot(To_Use_cells_Mar, main = "Cells to keep 'Marine'") #visualize distribution of 'valid' cells that remain 

##Extract ID's of usable/'valid' cells ----
Cells_to_use_Mar <- which(values(To_Use_cells_Mar)==T) #converts logical raster into a list of cell indices (numeric IDs) 

##Extract species composition for each grid cell ----
Phylacine_DataFrame_Marine_Poj_Binary <- values(Phylacine_rast_Marine_Poj_Binary) #converts raster stack into dataframe

## Spp_per_Use_cells_Mar ----
#Loop through each "usable/valid" grid cell
#Extract which species are present (value = 1)
#Match each species to its IUCN status from Trait_data_filtered
Spp_per_Use_cells_Mar <- lapply(Cells_to_use_Mar,
                            function(i){
                              Spp_List <- Phylacine_DataFrame_Marine_Poj_Binary[i,]
                              #Extract species presence/absence for this cell
                              Out <- data.frame(Spp=names(Spp_List[Spp_List!=0]))
                              #Keep only species that are present (value ≠ 0)
                              Out$IUCN <- Trait_data_filtered$IUCN.Status.1.2[match(Out$Spp,Trait_data_filtered$Binomial.1.2)]
                              #Match IUCN status using Trait_data_filtered and add as new column
                              return(Out)
                              #Return a dataframe: species name + IUCN category
                            })
names(Spp_per_Use_cells_Mar) <- paste0("Cell.ID.",Cells_to_use_Mar) #Name each list element by its corresponding cell ID
Spp_per_Use_cells_Mar[1:10] #Inspect the first 10 elements

#Output - Data frame (Spp_per_Use_cells_Mar)
  #List of data frames (cell ID's) 
    #For each cell ID a data frame: species name + matching IUCN category
      #$Spp - vector of species present in the grid
      #$IUCN - IUCN assessment of each species ('LC', 'VU', 'EN', 'CR')


#ANALYSIS Marine ----

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

#Spp_per_Use_cells_Mar = list of data frames, each element corresponds to one of the filtered grid cells (≤ 10 species + ≤ 2 EN species)
#Complete_trees = multiphylo object (list of 1000 full phylogenetic trees)


## Chunked tree Pre-Pruning ----
  #Needed do to very large amount of cells (+22.000 cells)
    #Strategy: Pre-prune trees in small batches and save to disk
      #only load what's needed during analysis

### Step 1: Create pruned tree Cache ----
  #Run once before analysis

#Creates pruning function
create_pruned_tree_cache <- function(cell_data, 
                                     tree_list,
                                     cache_dir = "pruned_trees_cache",
                                     batch_size = 100) {
  
  #Create cache directory
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir)
  }
  
  #Identify unique species compositions
    #Loop through all cells and create signatures from present species assemblages
    ##Create data frame and add column with signatures
  species_signatures <- sapply(cell_data, function(x) {
    paste(sort(x$Spp), collapse = "|")
  })
  
    #Sort through signatures and keep only unique 
  unique_signatures <- unique(species_signatures)
    #Amount of different signatures
  n_unique <- length(unique_signatures)
  
  #Print analysis parameters to console (Not needed - excess)
  cat(sprintf("Found %d unique species compositions\n", n_unique))
  cat(sprintf("Will create %d pruned tree sets (instead of %d)\n", 
              n_unique, length(cell_data)))
  cat("================================================================================\n") #create section divide to improve 'readability' when printed in console
  
  #Create mapping file (signature -> cell indices)
  mapping <- list()
  for (sig in unique_signatures) {
    mapping[[sig]] <- which(species_signatures == sig)
  }
  
  #Save mapping
  saveRDS(mapping, file.path(cache_dir, "signature_mapping.rds"))
  saveRDS(species_signatures, file.path(cache_dir, "species_signatures.rds"))
  
  #Calculation of batch size
    #Process in batches to avoid memory overload
  n_batches <- ceiling(n_unique / batch_size)
  
  #Print analysis update to console (Not needed - excess)
  cat(sprintf("\nProcessing %d unique compositions in %d batches of ~%d\n",
              n_unique, n_batches, batch_size))
  
  start_time <- Sys.time() #Track processing time
  
  #Loop through batchs
  for (batch_idx in 1:n_batches) {
    #Calculate which signatures are in this batch
    batch_start <- (batch_idx - 1) * batch_size + 1
    batch_end <- min(batch_idx * batch_size, n_unique)
    batch_sigs <- unique_signatures[batch_start:batch_end] #Creates vector of signature indices to process
    
    #Print analysis update to console (Not needed - excess)
    cat(sprintf("Batch %d/%d: Processing signatures %d-%d...\n",
                batch_idx, n_batches, batch_start, batch_end))
    
    batch_time_start <- Sys.time() #Track batch processing time
    
    #Process each signature in this batch
    for (sig_idx in seq_along(batch_sigs)) {
      sig <- batch_sigs[sig_idx]
      
      #Get species list for this signature
      cell_idx <- mapping[[sig]][1]  #Get first cell with this signature
      species_list <- cell_data[[cell_idx]]$Spp
      
      #Loop through and prune all trees for this species composition
      pruned_trees <- lapply(tree_list, function(tree) {
        species_to_keep <- species_list[species_list %in% tree$tip.label]
        
        if (length(species_to_keep) > 1 && 
            length(species_to_keep) < length(tree$tip.label)) {
          drop.tip(tree, tree$tip.label[!tree$tip.label %in% species_to_keep])
        } else {
          tree  #Return original if can't prune
        }
      })
      
      #Save this pruned tree set (as multiphylo object)
      class(pruned_trees) <- "multiPhylo"
      
      #Create safe filename from signature index
      sig_number <- batch_start + sig_idx - 1
      filename <- file.path(cache_dir, sprintf("pruned_trees_%05d.rds", sig_number))
      saveRDS(pruned_trees, filename)
      
      #Print progress within batch (Not needed - excess)
      if (sig_idx %% 10 == 0) {
        cat(sprintf("  Completed %d/%d in batch\n", sig_idx, length(batch_sigs)))
      }
    } #End of batch
    
    #Calculate batch processing time
    batch_time <- as.numeric(Sys.time() - batch_time_start, units = "secs")
    #Print process time for batch - 'progress bar' (Not needed - excess)
    cat(sprintf("Batch %d complete in %.1f seconds (%.2f sec/composition)\n\n",
                batch_idx, batch_time, batch_time / length(batch_sigs)))
    
    #Force garbage collection between batches
    gc(verbose = FALSE)
  }
  
  #Calculate total processing time
  total_time <- as.numeric(Sys.time() - start_time, units = "mins")
  
  #Display completion summary
  cat("================================================================================\n") #create section divide to improve 'readability' when printed in console
  cat(sprintf("Total time: %.2f minutes (%.2f hours)\n", total_time, total_time/60)) #time in minutes
  cat(sprintf("Cache location: %s/\n", cache_dir)) #File location (used for the analysis)
  cat(sprintf("Unique compositions: %d\n", n_unique)) #Number of unique species compositions
  
  #Calculate cache size
  cache_files <- list.files(cache_dir, pattern = "pruned_trees_.*\\.rds", full.names = TRUE)
  total_size_mb <- sum(file.info(cache_files)$size) / 1024^2
  cat(sprintf("Total cache size: %.1f MB\n", total_size_mb)) #Size of file
  
  #Returns cache directory without printing it in the console
  return(invisible(cache_dir))
}


### Step 2: Analysis function using pruned tree cache ----
  #Creates analysis function
calculate_ses_pd_cached <- function(cell_data,
                                    cache_dir = "pruned_trees_cache",
                                    n_rand = 1000,
                                    output_file = "SES_Summ_Final.csv",
                                    n_cores = NULL,
                                    chunk_size = 50) {
  
  #Pakages needed for analysis - in basic script
  library(parallel)
  library(foreach)
  library(doParallel)
  library(picante)
  library(ape)
  
  #Setup parallel backend
  if (is.null(n_cores)) {
    #If number of cores aren't specified it detects number of available cores minus 1 (keeping computer responsive)
    n_cores <- max(1, detectCores(logical = FALSE)) #leave out logical cores
  }
  
  #Print analysis parameters to console (Not needed - excess)
  cat("=== Caches tree analysis ===\n")
  cat(sprintf("Cells: %d | Cores: %d | Randomizations: %d\n",
              length(cell_data), n_cores, n_rand))
  cat("================================================================================\n") #create section divide to improve 'readability' when printed in console
  
  #Print analysis update to console (Not needed - excess)
  cat("Loading cache metadata...\n")
  #Load mapping (pre prunned cache)
  mapping <- readRDS(file.path(cache_dir, "signature_mapping.rds"))
  species_signatures <- readRDS(file.path(cache_dir, "species_signatures.rds"))
  
  #Create reverse lookup: cell index -> signature number
    #Assigns relevant signature to cell (which to use for analysis)
  sig_to_number <- list()
  sig_list <- names(mapping)
  for (i in seq_along(sig_list)) {
    sig_to_number[[sig_list[i]]] <- i
  }
  
  #Print analysis update to console (Not needed - excess)
  cat("Cache loaded. Starting analysis...\n\n")
  
  #Create cluster of worker processes (one pr. core)
  cl <- makeCluster(n_cores)
  #Register cluster so foreach knows to use it
  registerDoParallel(cl)
  
  #Export required packages to workers
  clusterEvalQ(cl, {
    library(picante) #Needed for pd() function
    library(ape)     #Needed for drop.tip() function (pruning)
  })
  
  #Export cache info to workers
  clusterExport(cl, c("cache_dir", "species_signatures", "sig_to_number"),
                envir = environment())
  
  #Initialize output file
    #If file exists from previous run delete it and start fresh
  if (file.exists(output_file)) file.remove(output_file)
  
  #Process cells in chunks to manage memory
  n_cells <- length(cell_data) #Total number of grid cells (21.819 for marine)
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
                                 Comm <- cell_data[[i]] #Extract species assemblage for this cell
                                 
                                 #Calculate richness
                                 total_richness <- nrow(Comm) #Species richness in the cell
                                 en_cr_richness <- sum(Comm$IUCN %in% c("EN", "CR")) #Number of Endangered species in this cell
                                 
                                 #Get species names
                                 spp_names <- make.unique(as.character(Comm$Spp))
                                 
                                 #Pre-allocate random community matrix
                                 n_species <- total_richness #Number of species in cell
                                 n_to_remove <- en_cr_richness #Number of species to remove in each random community
                                 #Create Matrix: rows = species, columns = random communities (1-1000), all values = 1
                                 rand_matrix <- matrix(1, nrow = n_species, ncol = n_rand) #n_rand - number of random communities wanted
                                 
                                 #Vectorized random community generation
                                  #For each random community, remove the same number of species, as number of EN/CR present, but randomly
                                 if (n_to_remove > 0 && n_to_remove < n_species) {
                                   rand_indices <- replicate(n_rand, sample.int(n_species, n_to_remove))
                                   for (j in 1:n_rand) { #Loop through each random community
                                     rand_matrix[rand_indices[, j], j] <- 0 #Set choosen species value to 0 (absent) in this random community
                                   }
                                 }
                                 
                                 #Build community matrix
                                  #Create presence/absence columns for different scenarios (original species assemblage, species assemblage without endangered)
                                  #Combine original and random communities
                                 comm_matrix <- cbind(
                                   Org = rep(1, n_species),
                                   No.EN = ifelse(Comm$IUCN %in% c("EN", "CR"), 0, 1),
                                   rand_matrix
                                 )
                                 rownames(comm_matrix) <- spp_names
                                 comm_matrix_t <- t(comm_matrix)
                                 #Columns: Org (all present), No.EN (EN/CR removed), Rand.1, Rand.2, ..., Rand.1000
                                 rownames(comm_matrix_t) <- make.unique(c("Org", "No.EN", 
                                                                          paste0("Rand.", 1:n_rand)))
                                 
                                 #Load pre-pruned trees for this cell
                                 sig <- species_signatures[i]
                                 sig_number <- sig_to_number[[sig]]
                                 pruned_file <- file.path(cache_dir, sprintf("pruned_trees_%05d.rds", sig_number))
                                 tree_list <- readRDS(pruned_file)
                                 
                                 #Calculate PD for all trees
                                  #Pre-allocate vectors to store results from each tree
                                 n_trees <- length(tree_list) #Number of trees in multiphylo
                                  #Create empty numeric vectors
                                 ses_values <- numeric(n_trees)
                                 pd_org_all <- numeric(n_trees)
                                 pd_no_en_all <- numeric(n_trees)
                                 mean_pd_rnd_all <- numeric(n_trees)
                                 sd_pd_rnd_all <- numeric(n_trees)
                                 
                                 #Calculate PD (trees already pruned)
                                  #Loop through each phylogenetic tree (d is the tree index)
                                 for (d in 1:n_trees) {
                                   #calculate PD for each community in matrix
                                    #include.root = FALSE - measures from MRCA of present species, not tree root
                                   pd_vals <- pd(comm_matrix_t, tree_list[[d]], include.root = FALSE)
                                   
                                   #Extract values
                                    #pd_vals dataframe contain rows for each community and columns $PD and $SR
                                   pd_org <- pd_vals$PD[1] #PD of row 1 (Original community)
                                   pd_no_en <- pd_vals$PD[2] #PD of row 2 (No EN/CR community)
                                   pd_rand <- pd_vals$PD[3:nrow(pd_vals)] #PD of rows 3 onward (all random communities)pd_rand - vector of PD results)
                                   
                                   #Calculate losses
                                    #How much PD is lost when species are removed
                                   observed_loss <- pd_org - pd_no_en #loss from removing EN/CR
                                   random_losses <- pd_org - pd_rand  #losses from removing random species (vector)
                                   
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
                                 data.frame(
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
                                 
                               }, #End of tryCatch block
                               error = function(e) { #(Extra safety)
                                 #If anything goes wrong in the try block above, this code runs instead
                                  #Returns error information for debugging
                                 data.frame(
                                   Cell.ID = names(cell_data)[i],
                                   Species.Richness = NA, EN.CR.Richness = NA,
                                   SES.Mean = NA, SES.SD = NA,
                                   PD.Org.Mean = NA, PD.Org.SD = NA,
                                   PD.NoEN.Mean = NA, PD.NoEN.SD = NA,
                                   PD.Rnd.Mean = NA, PD.Rnd.SD = NA,
                                   #stores error message
                                   Error = paste("ERROR:", as.character(e)),
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
    elapsed_total <- as.numeric(Sys.time() - overall_start, units = "secs") #Time (seconds) it took to analyse 'this far'
    cells_processed <- chunk_end #Total number of cells completed so far
    cells_remaining <- n_cells - cells_processed #Number of cells left to process
    avg_per_cell <- elapsed_total / cells_processed #Average seconds per cell in this chunk
    eta <- cells_remaining * avg_per_cell #Rough estimate of remaining time
    
    #Print progress bar to console
    cat(sprintf("Chunk %d/%d | Cells: %d-%d/%d (%.1f%%) | %.1fs | Avg: %.2fs/cell | ETA: %.1f min (%.1f hrs)\n",
                chunk_idx, n_chunks, chunk_start, chunk_end, n_cells,
                (cells_processed/n_cells)*100, chunk_time,
                avg_per_cell, eta/60, eta/3600)) #estimate time left in minutes and hours
    
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
  
  #Return csv file
  return(read.csv(output_file))
} #End of function


### Function usage - workflow ----

#### Step 1: Pre-pruning cache ----
cache_dir <- create_pruned_tree_cache(
  cell_data = Spp_per_Use_cells_Mar,         #Input: List of Marine cells
  tree_list = Complete_trees,                #Input: 1000 phylogenetic trees
  cache_dir = "marine_pruned_trees_cache",   #filename - save results to this file
  batch_size = 100                           #process 100 unique compositions at a time
)


#### Step 2: Run analysis using cache ----
results_marine_full <- calculate_ses_pd_cached(
  cell_data = Spp_per_Use_cells_Mar,         #Input: List of Marine cells
  cache_dir = "marine_pruned_trees_cache",   #Input: Cache of pre prunned trees
  n_rand = 1000,                             #1000 random communities per cell
  output_file = "SES_Summ_Final_Mar.csv",    #filename - save results to this file
  n_cores = detectCores(logical = FALSE),    #Number of cores used for analysis
  chunk_size = 50                            #Number of cells in chunks
)

#Convert csv file into Data frame
  #results_mar <- read.csv("SES_Summ_Final_Mar.csv")
if(!"Results_mar_df.rds" %in% dir()){
  saveRDS(results_mar, file = "Results_mar_df.rds") #saves data frame as rds file (faster to load)
} else {
  results_mar <- readRDS("Results_mar_df.rds") #Load rds file if file already exist
}


### Output file structure ----
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
#FINAL RUN
# Full data (1000 trees and 1000 Random communities using 8 cores)
# Output name:  "SES_Summ_Final_Mar.csv"
# Runtime:       145.9 hours
# Dates:         30.11.25 - 11.15 am
#                07.12.25 - 17.00 pm
