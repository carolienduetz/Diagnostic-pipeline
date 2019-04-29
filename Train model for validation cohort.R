source("helperFunctions_diagnostic_pipeline.R")

# Create final model -----------------------------------------------------------

outputDir <- "Final_model"
dir.create(outputDir)
version <-""

# parameters -------------------------------------------------------------------
verbose <- FALSE
seeds <- 1
seed <- 1
nFolds <- 5
percentTestPerFold <- 0.2
tubes <- 2:7
patient_labels <- c("1" = "MDS",
                    "2" = "Control")

# FlowSOM parameters
cols_to_use <- c(1, 3:7, 9:12)
nCellsPerFile <- 40000
gridSize <- 15
nClus <- 30
nFeatures_options <- ("All")

# Classification parameters
feature_types <- c(#pC= "^pC[0-9]", 
                   #rpC = "^rpC[0-9]*", 
                   MC = "^MC[0-9]*")

FS_method <- "randomForest" # "randomForest" or "mRMR"
classification_method <- "randomForest" #GMnet" # "randomForest" or "SVM"

# read metadata

metadata <- read_metadata("") %>% 
  dplyr::mutate(MDS_class = as.factor(patient_labels[as.character(MDS)])) %>% 
  dplyr::filter(!is.na(MDS_class))

# build FlowSOM    -------------------------------------------------------------

for(tube in tubes){
  preprocessed_dir <- paste0("preprocessed tube ", tube)
  files_in_dir <- get_files_in_dir(preprocessed_dir)
  
  files_final <- files_in_dir[ metadata %>% 
                                 dplyr::filter(Fenotypering_nr %in% names(files_in_dir)) %>% 
                                 dplyr::pull(Fenotypering_nr)]
  
  fsommodel <- FlowSOM_training(metadata = metadata,
                            files_train = files_final,
                            cols_to_use = cols_to_use,
                            nCellsPerFile = nCellsPerFile,
                            gridSize = gridSize,
                            nClus = nClus,
                            seed = seed,
                            verbose = TRUE)

  
  model_file <- file.path(outputDir,
                          paste0("final_fsommodel_",version, tube,".RDS"))
  saveRDS(fsommodel, 
          file = model_file)
}

# extract features from FlowSOM tree -------------------------------------------

full_data_new <- NULL
for(tube in tubes){
  
  
  model_file <- file.path(outputDir,
                          paste0("final_fsommodel_",version, tube,".RDS"))
  fsommodel <- readRDS(model_file)
  
  features_to_use <- FS_byName(fsommodel$full_data, feature_types)
  feature_ids <- which(colnames(fsommodel$full_data) %in% features_to_use)
  colnames(fsommodel$full_data)[feature_ids] <-  
    paste0("T",
           tube,
           "_",
           colnames(fsommodel$full_data)[feature_ids])
  if (is.null(full_data_new)) {
    full_data_new <- fsommodel$full_data
    
  } else {
    full_data_new <- full_join(full_data_new, 
                                 fsommodel$full_data[,c("Fenotypering_nr",
                                                    paste0("T", tube, "_", features_to_use))], 
                                 by = "Fenotypering_nr")
  }
}


rownames(full_data_new) <- full_data_new$Fenotypering_nr


# Feature selection ------------------------------------------------------------
  
features_to_use <- FS_byName(full_data_new,
                             gsub("\\^", "^T[0-9]*_", feature_types))

  
  if(FS_method == "mRMR"){
    # MRMR preparation to select multiple feature subsets
    dd <- mRMR.data(data = cbind(as.numeric(full_data_train[, "MDS_class"]),
                                 full_data_train[, features_to_use]))
  } else if(FS_method == "randomForest"){
    # RF preparation to select multiple feature subsets
    forest <- RF_training(data = full_data_new,
                          features_to_use = features_to_use,
                          num.trees = 50000,
                          seed = seed)
  }
  
  for(nFeatures in nFeatures_options){
    if(nFeatures == "All"){
      selected_features <- features_to_use
    } else {
      
      if(FS_method == "mRMR"){
        # MRMR feature selection
        filter <- mRMR.classic(data = dd,
                               target_indices = 1,
                               feature_count = as.numeric(nFeatures))
        selected_features <- features_to_use[rev(filter@filters[[1]]) - 1]
      } else if (FS_method == "randomForest") { 
        # Random forest feature selection
        selected_features <- FS_randomforest(forest, nFeatures)
      } else {
        stop("Only randomForest or mRMR are currently supported for feature ",
             "selection.")
      }
      
    }
    
# train classifier -------------------------------------------------------------
    
    prediction_file <- file.path(outputDir,
                                 paste0("prediction_", hoi,
                                        "_nF", nFeatures, ".RDS"))
    if (classification_method == "randomForest"){
      forest <- RF_training(data = full_data_new,
                            features_to_use = selected_features,
                            num.trees = 10000,
                            seed = seed)
      model_classifier <- list(features_to_use = selected_features,
                               classifier = forest)
    } else if (classification_method == "SVM") {
      svm <- SVM_training(data = full_data_new,
                          features_to_use = selected_features,
                          seed = seed)
      model_classifier <- list(features_to_use = selected_features,
                               classifier = svm)
    } else if (classification_method == "GMnet") {
      gmnet <- GMnet_training(data = full_data_train,
                              features_to_use = selected_features,
                              seed = seed)
      model_classifier <- list(features_to_use = selected_features,
                               classifier = gmnet)
    } else {
      stop("Only randomForest or SVM are currently supported for classification")
    }
    
  }
    final_model <- model_classifier
    
    saveRDS(final_model,
            paste0(""))
    