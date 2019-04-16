
# load functions from other file -----------------------------------------------

source("helperFunctions_diagnostic_pipeline.R")

# Set some parameters ----------------------------------------------------------

outputDir <- "results_MDSvsControl/"
dir.create(outputDir)

verbose <- FALSE
seeds <- 3
nFolds <- 5
percentTestPerFold <- 0.2
tubes <- 4
patient_labels <- c(
  "1" = "MDS",
  "2" = "Control")

# FlowSOM parameters -----------------------------------------------------------
cols_to_use <- c(1, 3:7, 9:12)
nCellsPerFile <- 40000
gridSize <- 15
nClus <- 30
nFeatures_options <- c("All")

# Classification parameters
feature_types <- c(
  pC= "^pC[0-9]*", 
  rpC = "^rpC[0-9]*", 
  MC = "^MC[0-9]*")

FS_method <- "randomForest" #"mRMR" # or  
classification_method <- "randomForest" # "GMnet" #""SVM"  # """"#"" or  

parameter_string <- paste0("Markers_",paste(cols_to_use, collapse="_"),
                           "_cells", nCellsPerFile,
                           "_grid", gridSize,
                           "_clus", nClus,
                           "_features", paste0(names(feature_types), collapse="_"))

# Read the metadata ------------------------------------------------------------


metadata <- read_metadata("") %>% 
  dplyr::mutate(MDS_class = as.factor(patient_labels[as.character(MDS)])) %>% 
  dplyr::filter(!is.na(MDS_class))


#dplyr::filter(T10 > nCells)

# Split in folds ---------------------------------------------------------------

metadata$Fold <- as.factor(as.numeric(cut(as.numeric(as.factor(metadata$Fenotypering_nr)), 
                                          nFolds)))
ggplot(metadata) +
  geom_point(aes(x = MDS, y = Fenotypering_nr, col = Fold)) +
  theme_minimal()


# Compute for every tube seperately --------------------------------------------


aucs <- matrix(0,
               nrow = length(seeds) * nFolds,
               ncol = 1+length(tubes),
               dimnames = list(apply(expand.grid(seeds, seq_len(nFolds)), 1, 
                                     paste, collapse="_"),
                               c(paste0("Tube ",tubes), "All tubes")))
aucs <- rep(list(aucs), length(nFeatures_options))
names(aucs) <- as.character(nFeatures_options)


seed <- seeds[1]

for(seed in seeds){
  fold <- nFolds[1]
  for(fold in seq_len(nFolds)){
    tube <- tubes[1]
    for(tube in tubes){
      
      message("Processing tube ",tube)
      settings_string <- paste0("Tube",tube,
                                "_seed", seed,
                                "_fold", fold,
                                parameter_string)
      
      preprocessed_dir <- paste0("preprocessed tube ", tube)
      files_in_dir <- get_files_in_dir(preprocessed_dir)  
      
      files_train <- files_in_dir[ metadata %>% 
                                     dplyr::filter(Fenotypering_nr %in% names(files_in_dir)) %>% 
                                     dplyr::filter(Fold != fold) %>%    
                                     dplyr::pull(Fenotypering_nr)]
      
      files_test <- files_in_dir[ metadata %>% 
                                    dplyr::filter(Fenotypering_nr %in% names(files_in_dir)) %>% 
                                    dplyr::filter(Fold == fold) %>%      
                                    dplyr::pull(Fenotypering_nr)]
      
      model_file <- file.path(outputDir,
                              paste0("model_",settings_string,".RDS"))
      
      # FlowSOM ----------------------------------------------------------------
      
      if (file.exists(model_file)) {
        model <- readRDS(model_file)
        message("  Reloaded FlowSOM model.")
      } else {
        model <- FlowSOM_training(metadata = metadata,
                                  files_train = files_train,
                                  cols_to_use = cols_to_use,
                                  nCellsPerFile = nCellsPerFile,
                                  gridSize = gridSize,
                                  nClus = nClus,
                                  seed = seed,
                                  verbose = TRUE)
        saveRDS(model, 
                file = model_file)
      }
      
      PlotStars(model$FlowSOM$FlowSOM)
      
      # Feature selection ------------------------------------------------------
      
      features_to_use <- FS_byName(model$full_data,
                                   feature_types)
      if(FS_method == "mRMR"){
        # MRMR preparation to select multiple feature subsets
        dd <- mRMR.data(data = cbind(as.numeric(model$full_data[, "MDS_class"]),
                                     model$full_data[, features_to_use]))
      } else if(FS_method == "randomForest"){
        # RF preparation to select multiple feature subsets
        forest <- RF_training(data = model$full_data,
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
        
        # Train classifier -----------------------------------------------------
        
        if (classification_method == "randomForest"){
          forest <- RF_training(data = model$full_data,
                                features_to_use = selected_features,
                                num.trees = 10000,
                                seed = seed)
          model_classifier <- list(features_to_use = selected_features,
                                   classifier = forest)
        } else if (classification_method == "SVM") {
          svm <- SVM_training(data = model$full_data,
                              features_to_use = selected_features,
                              seed = seed)
          model_classifier <- list(features_to_use = selected_features,
                                   classifier = svm)
        } else if (classification_method == "GMnet") {
          gmnet <- GMnet_training(data = model$full_data,
                                  features_to_use = selected_features,
                                  seed = seed)
          model_classifier <- list(features_to_use = selected_features,
                                   classifier = gmnet)
        } else {
          stop("Only randomForest or SVM are currently supported for classification")
        }
        
        # Test dataset on classifier -------------------------------------------
        
        prediction_file <- file.path(outputDir,
                                     paste0("prediction_", settings_string,
                                            "_nF", nFeatures, ".RDS"))
        
        if (file.exists(prediction_file)) {
          prediction <- readRDS(prediction_file)
          message("  Reloaded prediction.")
        } else {
          
          new_data <- new_files_to_FlowSOM(model$FlowSOM,
                                           files_test,
                                           metadata,
                                           verbose = TRUE)
          
          prediction <- model_testing(model = model_classifier,
                                      data = new_data,
                                      classification_method = classification_method)
          
          saveRDS(prediction, file = prediction_file)
        }
        
        predicted_labels <- c("Control", "MDS")[1 + (prediction$prediction[,2] >= 0.5)] 
        
        real_labels <- metadata %>% 
          dplyr::filter(Fenotypering_nr %in% names(files_test)) %>%
          dplyr::pull("MDS_class")
        
        confusion_matrix <- table(real_labels, predicted_labels)
        print(confusion_matrix)
        
        toWrite <- rbind(confusion_matrix,
                         "",
                         c(paste(prediction$full_data$Fenotypering_nr[which(predicted_labels != real_labels)],
                                 collapse = ", "),
                           ""))
        if(nFeatures == "All") { 
          write.csv2(toWrite,
                     paste0("confusion_",settings_string,"_nF",nFeatures,".csv"))
        }
        
        message("Wrongly predicted patients: ",
                paste(prediction$data$Fenotypering_nr[which(predicted_labels != real_labels)], collapse = ", "))
        
        aucs[[as.character(nFeatures)]][paste0(seed, "_", fold), 
                                        paste0("Tube ",tube)] <- performance(prediction$evaluation, "auc")@y.values[[1]]
        
      }
      
    } 
    
    
    # Combine all tubes --------------------------------------------------------
    full_data_train <- NULL
    full_data_test <- NULL
    
    for(tube in tubes){
      
      settings_string <- paste0("Tube",tube,
                                "_seed", seed,
                                "_fold", fold,
                                parameter_string)
      
      model_file <- file.path(outputDir,
                              paste0("model_", settings_string,".RDS"))
      model <- readRDS(model_file)
      
      prediction_file <- file.path(outputDir,
                                   paste0("prediction_", settings_string,"_nF",nFeatures,".RDS"))
      prediction <- readRDS(prediction_file)
      
      features_to_use <- FS_byName(model$full_data, feature_types)
      feature_ids <- which(colnames(model$full_data) %in% features_to_use)
      colnames(model$full_data)[feature_ids] <-  
        colnames(prediction$data)[feature_ids] <- 
        paste0("T",
               tube,
               "_",
               colnames(model$full_data)[feature_ids])
      
      if (is.null(full_data_train)) {
        full_data_train <- model$full_data
        full_data_test <- prediction$data
      } else {
        full_data_train <- full_join(full_data_train, 
                                     model$full_data[,c("Fenotypering_nr",
                                                        paste0("T", tube, "_", features_to_use))], 
                                     by = "Fenotypering_nr")
        full_data_test <- full_join(full_data_test, 
                                    prediction$data[,c("Fenotypering_nr",
                                                       paste0("T", tube, "_", features_to_use))], 
                                    by = "Fenotypering_nr")
      }
    }
    rownames(full_data_train) <- full_data_train$Fenotypering_nr
    rownames(full_data_test) <- full_data_test$Fenotypering_nr
    
    saveRDS(list(full_data_train, full_data_test),
            file = file.path(outputDir,
                             paste0("dataAllTubes_",settings_string,".RDS")))
    
    # Prepare data for feature selection ---------------------------------------
    features_to_use <- FS_byName(full_data_train,
                                 gsub("\\^", "^T[0-9]*_", feature_types))
    
    if(FS_method == "mRMR"){
      # MRMR preparation to select multiple feature subsets
      dd <- mRMR.data(data = cbind(as.numeric(full_data_train[, "MDS_class"]),
                                   full_data_train[, features_to_use]))
    } else if(FS_method == "randomForest"){
      # RF preparation to select multiple feature subsets
      forest <- RF_training(full_data_train,
                            features_to_use,
                            num.trees = 50000,
                            seed = seed)
    }
    
    # Perform feature selection ------------------------------------------------
    for(nFeatures in nFeatures_options){
      
      if(nFeatures == "All"){
        selected_features <- features_to_use
      } else { 
        
        if(FS_method == "mRMR"){
          # MRMR feature selection ---------------------------------------------
          filter <- mRMR.classic(data = dd,
                                 target_indices = 1,
                                 feature_count = as.numeric(nFeatures))
          selected_features <- features_to_use[rev(filter@filters[[1]]) - 1]
          
        } else if (FS_method == "randomForest") { 
          
          # Random forest feature selection ------------------------------------
          selected_features <- FS_randomforest(forest, as.numeric(nFeatures))
          
        } else {
          stop("Only randomForest or mRMR are currently supported for feature ",
               "selection.")
        }
        
        # Save selected features to RDS ----------------------------------------
        saveRDS(selected_features,
                file = file.path(outputDir,
                                 paste0("selectedAllTubes_",settings_string,"_nF",nFeatures,".RDS")))
      }
      
      # Train classification model --------------------------------------------- 
      
      if (classification_method == "randomForest"){
        forest <- RF_training(data = full_data_train,
                              features_to_use = selected_features,
                              num.trees = 10000,
                              seed = seed)
        model_classifier <- list(features_to_use = selected_features,
                                 classifier = forest)
      } else if (classification_method == "SVM") {
        svm <- SVM_training(data = full_data_train,
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
      
      # Predict ----------------------------------------------------------------
      
      prediction <- model_testing(model = model_classifier,
                                  data = full_data_test,
                                  classification_method = classification_method)
      
      saveRDS(prediction, file= file.path(outputDir, 
                                          paste0("predictionAllTubes_", settings_string, ".RDS")))
      
      COMP <- prediction(prediction$prediction[,2], full_data_test[,"MDS_class"])
      predicted_labels <- c("Control", "MDS")[1 + (prediction$prediction[,2] >= 0.5)]
      real_labels <- metadata %>% 
        dplyr::filter(Fenotypering_nr %in% names(files_test)) %>%
        dplyr::pull("MDS_class")
      table(real_labels, predicted_labels)
      
      confusion_matrix <- table(real_labels, predicted_labels)
      print(confusion_matrix)
      
      message("Wrongly predicted patients: ",
              paste(prediction$data$Fenotypering_nr[which(predicted_labels != real_labels)], collapse = ", "))
      
      toWrite <- rbind(confusion_matrix,
                       "",
                       c(paste(names(which(predicted_labels != real_labels)), collapse = ", "),""))
      if(nFeatures == "All") { 
        write.csv2(toWrite,
                   paste0("confusion_allTubes_",settings_string,"_nF",nFeatures,".csv"))
      }
      
      aucs[[as.character(nFeatures)]][paste0(seed, "_", fold), 
                                      "All tubes"] <- performance(COMP, "auc")@y.values[[1]]
      
    }
  }
}

boxplot(aucs, main = paste0("AUC ",nFolds,"-fold cross validation"))

