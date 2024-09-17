library(vegan)
library(dplyr)
library(phytools)
library(phangorn)

Grass_Features <- read.csv("~/Desktop/Grass_Images_Patches_Modern_Concatenated_Features.csv", header=FALSE)
Grass_Labels_All <- read.csv("~/Desktop/Modern_Poaceae_Specimens_Labels.csv", header=FALSE)

Labels <- Grass_Labels_All$V1
Numeric_Labels =  Labels + 1

############################

Grass_Features$Label <- Numeric_Labels # Last column now contains labels (species ID)

# Step 2 & 3: Group by the label and calculate the mean for each group
average_features <- Grass_Features %>%
  group_by(Label) %>%
  summarise(across(everything(), ~mean(.x, na.rm = TRUE)))

average_grass_features <- average_features[2:4097]


species_mapping <- c('Agropyron_ciliare_Elymus_ciliaris', 'Agrostis_mexicana_Agrostis_tolucensis', 'Agrostis_quinqueseta', 'Agrostis_trachyphylla', 'Agrostis_volkensii', 'Amphibromus_neesii', 'Andropogon_amethystinus', 'Andropogon_chrysostachyus', 'Andropogon_contortus_Heteropogon_contortus', 'Andropogon_lima', 'Andropogon_schirensis', 'Anthoxanthum_nivale', 'Aristida_implexa_Aristida_megapotamica', 'Bothriochloa_bladhii', 'Brachiaria_brizantha_Urochloa_brizantha', 'Brachypodium_flexum', 'Bromus_auleticus', 'Bromus_brachyphyllus_Bromus_orcuttianus', 'Bromus_ciliatus', 'Bromus_lanatus', 'Bromus_leptoclados', 'Bromus_thominei_Bromus_hordeaceus_subsp_thominei', 'Calamagrostis_epigejos', 'Chrysopogon_fallax', 'Cymbopogon_nardus', 'Cynodon_dactylon', 'Dactylis_glomerata', 'Digitaria_abyssinica', 'Echinopogon_caespitosus', 'Ehrharta_erecta', 'Ehrharta_erecta_var_natalensis', 'Eleusine_coracana', 'Eleusine_jaegeri', 'Exotheca_abyssinica', 'Festuca africana_Pseudobromus_africanus', 'Festuca_costata', 'Festuca_elatior_Ampelodesmos_mauritanicus', 'Isachne_mauritiana', 'Koeleria_capensis', 'Lolium_perenne', 'Melica_onoei', 'Miscanthus_violaceus_Miscanthidium_violaceum', 'Oplismenus_compositus', 'Oplismenus_hirtellus', 'Pennisetum_longistylum_Cenchrus_clandestinus', 'Pentaschistis_borussica_Pentameris_borussica', 'Phalaris_arundinacea', 'Poa_anceps', 'Poa_leptoclada', 'Poa_schimperiana', 'Saccharum_arundinaceum_Tripidium_arundinaceum', 'Secale_cereale', 'Setaria_megaphylla', 'Sinarundinaria_alpina_Oldeania_alpina', 'Sorghum_halepense', 'Spartina_pectinata_Sporobolus_michauxianus', 'Stipa_compacta_Austrostipa_flavescens', 'Streblochaete_longiarista_Koordersiochloa_longiarista', 'Themeda_triandra', 'Zea_mays')

rownames(average_grass_features) = species_mapping
averaged_grass_features_standardized <- decostand(average_grass_features, "standardize")

## Replace off-diagonal elements with epsilon

replace_off_diagonal_zeros_with_constant <- function(matrix) {
  # Get the dimensions of the matrix
  n <- nrow(matrix)
  # Identify the diagonal elements
  diag_indices <- diag(n) == 1
  # Replace off-diagonal zeros with 1
  matrix[matrix == 0 & !diag_indices] <- 1e-8
  # Return the modified matrix
  return(matrix)
}

Poaceae_distance_matrix <- replace_off_diagonal_zeros_with_constant(Poaceae_distance_matrix)
dist = as.dist(Poaceae_distance_matrix)

tree = upgma(dist, "average")
Poaceae_tree = as.phylo(tree)
pollen_pca <- phyl.pca(Poaceae_tree, average_grass_features, method="BM", mode="corr")

# Retrieve mean and standard deviation of modern Poaceae pollen data

modern_mean <- apply(average_grass_features, 2, mean)
modern_std <- apply(average_grass_features, 2, sd)

modern_standardized <- sweep(average_grass_features, 2, modern_mean, "-")
modern_standardized <- sweep(modern_standardized, 2, modern_std, "/")
modern_standardized <- as.matrix(modern_standardized)

# Get modern pPCA scores
modern_pca_scores <- modern_standardized %*% pollen_pca$Evec[, 1:59]

######
######
######

# Initialize an empty list to store PCA scores for each component
modern_pca_scores_list <- list()

# Define the maximum number of components
max_components <- 59

# Loop through each component
for (i in 1:max_components) {
  # Compute PCA scores for the current component
  modern_pca_scores <- modern_standardized %*% pollen_pca$Evec[, 1:i]
  # Store the scores in the list
  modern_pca_scores_list[[paste("pPC", i, "_Modern", sep = "")]] <- modern_pca_scores[, i]
}

modern_pca_scores_list

max_components <- length(modern_pca_scores_list)
pca_columns <- setNames(modern_pca_scores_list, paste0('pPC', 1:max_components))

photosynthetic_pathway <- c(
  "C3",  # Agropyron_ciliare_Elymus_ciliaris
  "C4",  # Agrostis_mexicana_Agrostis_tolucensis_Changed to eragrostis and moved to C4
  "C3",  # Agrostis_quinqueseta
  "C3",  # Agrostis_trachyphylla
  "C3",  # Agrostis_volkensii
  "C3",  # Amphibromus_neesii
  "C4",  # Andropogon_amethystinus ..
  "C4",  # Andropogon_chrysostachyus
  "C4",  # Andropogon_contortus_Heteropogon_contortus
  "C4",  # Andropogon_lima
  "C4",  # Andropogon_schirensis_Changed to C4 _ Potential error in Wooller et al. 2001
  "C3",  # Anthoxanthum_nivale
  "C4",  # Aristida_implexa_Aristida_megapotamica
  "C4",  # Bothriochloa_bladhii
  "C4",  # Brachiaria_brizantha_Urochloa_brizantha
  "C3",  # Brachypodium_flexum
  "C3",  # Bromus_auleticus
  "C3",  # Bromus_brachyphyllus_Bromus_orcuttianus
  "C3",  # Bromus_ciliatus
  "C3",  # Bromus_lanatus
  "C3",  # Bromus_leptoclados
  "C3",  # Bromus_thominei_Bromus_hordeaceus_subsp_thominei
  "C3",  # Calamagrostis_epigejos
  "C4",  # Chrysopogon_fallax
  "C4",  # Cymbopogon_nardus
  "C4",  # Cynodon_dactylon
  "C3",  # Dactylis_glomerata
  "C4",  # Digitaria_abyssinica
  "C3",  # Echinopogon_caespitosus
  "C3",  # Ehrharta_erecta
  "C3",  # Ehrharta_erecta_var_natalensis _ Changed to C3..less likely to be C4
  "C4",  # Eleusine_coracana
  "C4",  # Eleusine_jaegeri
  "C4",  # Exotheca_abyssinica
  "C3",  # Festuca africana_Pseudobromus_africanus
  "C3",  # Festuca_costata
  "C3",  # Festuca_elatior_Ampelodesmos_mauritanicus
  "C3",  # Isachne_mauritiana
  "C3",  # Koeleria_capensis
  "C3",  # Lolium_perenne
  "C3",  # Melica_onoei
  "C4",  # Miscanthus_violaceus_Miscanthidium_violaceum
  "C3",  # Oplismenus_compositus
  "C3",  # Oplismenus_hirtellus
  "C4",  # Pennisetum_longistylum_Cenchrus_clandestinus
  "C3",  # Pentaschistis_borussica_Pentameris_borussica
  "C3",  # Phalaris_arundinacea
  "C3",  # Poa_anceps
  "C3",  # Poa_leptoclada
  "C3",  # Poa_schimperiana
  "C4",  # Saccharum_arundinaceum_Tripidium_arundinaceum
  "C3",  # Secale_cereale
  "C4",  # Setaria_megaphylla
  "C3",  # Sinarundinaria_alpina_Oldeania_alpina
  "C4",  # Sorghum_halepense
  "C4",  # Spartina_pectinata_Sporobolus_michauxianus
  "C3",  # Stipa_compacta_Austrostipa_flavescens
  "C3",  # Streblochaete_longiarista_Koordersiochloa_longiarista
  "C4",  # Themeda_triandra
  "C4"   # Zea_mays
)

combined_data <- data.frame(PhotosyntheticPathway = photosynthetic_pathway, pca_columns)

######
######
######

#### project fossil data onto pPCA (conducted using the modern data) and get their pPC scores 

# Load CNN feature matrix for the fossil data
Rutundu_Fossils <- read.csv("~/Desktop/Grass_Images_Patches_Fossils_Concatenated_Features.csv", header=FALSE)

# Standardize the fossil data using modern mean and std
Rutundu_Fossils_standardized <- sweep(Rutundu_Fossils, 2, modern_mean, "-")
Rutundu_Fossils_standardized <- sweep(Rutundu_Fossils_standardized, 2, modern_std, "/")
Rutundu_Fossils_standardized <- as.matrix(Rutundu_Fossils_standardized)

fossil_pca_scores <- Rutundu_Fossils_standardized %*% pollen_pca$Evec[, 1:59]

##################

library(caret)
library(randomForest)

photosynthetic_pathway <- as.factor(photosynthetic_pathway)


# Set up control for RFE with repeated cross-validation
rfeControl <- rfeControl(functions = rfFuncs,
                         method = "repeatedcv",
                         number = 5,
                         repeats = 10,
                         verbose = FALSE)

# Run RFE with Random Forest to find the optimal number of PCs
set.seed(123)  
rfeResults <- rfe(modern_pca_scores, 
                  photosynthetic_pathway, 
                  sizes = c(1:59),  # Test all 59 pPCs
                  rfeControl = rfeControl)

print(rfeResults)

# List of optimal PCs
optimal_pcs <- predictors(rfeResults)

# Train the final model using all modern species and the optimal pPCs (determined earlier using RFE)
finalModel <- randomForest(modern_pca_scores[, optimal_pcs], 
                           photosynthetic_pathway, 
                           ntree = 1000)

print(finalModel)

# Most important pPCs and their effects on the homogeneity of the nodes and leaves in the random forest
varImpPlot(finalModel, main = "Feature Importance")

# Model's accuracy vs. number of pPCs in random forest
plot(rfeResults$results$Variables, rfeResults$results$Accuracy, type = "b",
     xlab = "Variables", ylab = "Accuracy (Repeated Cross-Validation)",
     main = "")

# Select the same optimal pPCs from the fossil data
fossil_pca_scores_optimal <- fossil_pca_scores[, optimal_pcs]
# Predict photosynthetic pathway for each fossil specimen
fossil_predictions <- predict(finalModel, fossil_pca_scores_optimal, type = "prob")

print(fossil_predictions)

# Core depths
Y= c(20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,124,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,134,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,163,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,173,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,195,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,236,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,256,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,263,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,286,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,298,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,531,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,321,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,335,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,375,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,397,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,423,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,432,
     432,432,432,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,453,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,470,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,483,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516,516)

# Depth to Age mapping
depth_to_age <- c("20" = 252.922, "39" = 1006.79, "71" = 2446.333, "124" = 4497.479, 
                  "134" = 5001.556, "163" = 6496.673, "173" = 7022.491, "195" = 8195.662, 
                  "236" = 9281.874, "256" = 10554.34, "263" = 11004.19, "286" = 12497.5, 
                  "298" = 13285.34, "321" = 14230.91, "335" = 14995.26, "375" = 17204.36, 
                  "397" = 18434.3, "423" = 19520.83, "432" = 20040.09, "453" = 21257.91, 
                  "470" = 22249.94, "483" = 23012.15, "516" = 24014.17, "531" = 25000.81)

# Round ages to the nearest 5
rounded_depth_to_age <- round(depth_to_age / 5) * 5

ages <- sapply(Y, function(depth) rounded_depth_to_age[as.character(depth)])

ages_vector <- as.numeric(ages)  # Ensure it's numeric for any subsequent calculations

# Combine the predictions with their corresponding ages
fossil_data <- data.frame(Age = ages, Prediction = colnames(fossil_predictions)[max.col(fossil_predictions)],
                          C3_Prob = fossil_predictions[, "C3"], C4_Prob = fossil_predictions[, "C4"])

# Group by Age and calculate C4 relative abundance with bootstrapping
num_bootstraps <- 1000
bootstrap_results <- matrix(NA, nrow = num_bootstraps, ncol = length(unique(fossil_data$Age)))

# Perform bootstrapping
set.seed(123)
for (i in 1:num_bootstraps) {
  resample_indices <- sample(1:nrow(fossil_data), replace = TRUE)
  resample_data <- fossil_data[resample_indices, ]
  
  relative_abundance_resample <- aggregate(C4_Prob ~ Age, resample_data, mean)
  
  bootstrap_results[i, ] <- relative_abundance_resample$C4_Prob
}

# Calculate mean and confidence intervals
mean_abundance <- apply(bootstrap_results, 2, mean)
lower_ci <- apply(bootstrap_results, 2, function(x) quantile(x, 0.025))
upper_ci <- apply(bootstrap_results, 2, function(x) quantile(x, 0.975))

# Combine into a data frame
relative_abundance_with_ci <- data.frame(
  Age = unique(fossil_data$Age),
  Mean_C4_Relative_Abundance = mean_abundance,
  Lower_CI = lower_ci,
  Upper_CI = upper_ci
)

print(relative_abundance_with_ci)
