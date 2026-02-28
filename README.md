Project title:  
Effect of future extinctions of endangered species on phylogenetic diversity in terrestrial and marine mammals  

Description:  
Conatins all the raw data sourced from PHYLACINE version (1.2.1), the R scripts and the produced results, used in bachelor thesis.

Common Abriviations:  
terr - Terrestrial realm  
mar - Marine realm  
PD - Phylogenetic diversity  
SR - Species richness  
EN - Endangered species  

Overview: (fileformats)  
Analysis results  
- marine_pruned_trees_cache : Containes pre pruned multiphylos' stored by signature markers produced for marine analysis  
- SES_Summ_Final_Mar.csv : columns seperated by ',' contains results produced by Marine analysis  
- SES_Summ_Final_Terr.csv : columns seperated by ',' contains results produced by Terrestrial analysis

Phylacine 1.2.1  
- Current ranges : Contains range data for 5.832 mammal species (.tif)  
- Phylogenies : Contains multiphylo holding a 1000 possible phylogenies for all mammal species (.nex)
  (Uploaded as zip file due to GitHub upload size limitations)   
- Traits : Contains trait data (e.g. Taxonomy, Habitat, IUCN assessment...) (.csv - columns seperated by ',')  

Bassis_script.R (R script): Needed R packages, Data import and basis handling/filtering, vector creation (e.g. marine & terrestrial)  
Analysis_Mar.R (R script): Analysis of Marine Cells  
Analysis_Terr.R (R script): Analysis of Terrestrial Cells  
Figures.R (R script): Data vizualisation 

Requirements:  
Basis_script.R - Must be run before other scripts  
Figures.R - Requires that Basis_script + both Analysis scripts have been run  
