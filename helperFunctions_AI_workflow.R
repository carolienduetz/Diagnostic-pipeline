library(tidyverse)
library(magrittr)
library(flowCore)
library(FlowSOM)
library(randomForest)
library(rpart)
library(caret)
library(ROCR)
library(partykit)
library(mRMRe)
library(flowAI)

#' @param fsom  Result from calling FlowSOM function
#' @param files Named array with full path to training fcs files and 
#'              phenotype numbers as names
#' 
#' @return Data frame including counts (Cx), percentages (pCx), and 
#'         relative percentages (compared to the metacluster) per cluster (rCx) 
#'         and percentages per metacluster (MC)
new_files_to_FlowSOM <- function(fsom, 
                                 files, 
                                 metadata,
                                 MFI_markers = NULL,
                                 verbose = FALSE){
  
  if(is.null(MFI_markers)){ MFI_markers <- fsom$FlowSOM$map$colsUsed }
  
  # abbundance per node matrix ------------------------------------------------- 
  counts <- matrix(0,
                   nrow = length(files),
                   ncol = fsom$FlowSOM$map$nNodes,
                   dimnames = list(names(files),
                                   paste0("C",1:fsom$FlowSOM$map$nNodes)))
  MFI_clusters <- matrix(0,
                         nrow = length(files),
                         ncol = fsom$FlowSOM$map$nNodes * length(MFI_markers),
                   dimnames = list(names(files),
                                   apply(expand.grid(paste0("pC",1:fsom$FlowSOM$map$nNodes,"_MFI_"), 
                                                     fsom$FlowSOM$prettyColnames[MFI_markers]), 1 , paste, collapse="")))
  MFI_metaclusters <- matrix(0,
                             nrow = length(files),
                             ncol = max(as.numeric(fsom$metaclustering)) * length(MFI_markers),
                             dimnames = list(names(files),
                                             apply(expand.grid(paste0("MC",1:max(as.numeric(fsom$metaclustering)),"_MFI_"), 
                                                               fsom$FlowSOM$prettyColnames[MFI_markers]), 1 , paste, collapse="")))
  
  CV_clusters <- matrix(0,
                         nrow = length(files),
                         ncol = fsom$FlowSOM$map$nNodes * length(MFI_markers),
                         dimnames = list(names(files),
                                         apply(expand.grid(paste0("pC",1:fsom$FlowSOM$map$nNodes,"_CV_"), 
                                                           fsom$FlowSOM$prettyColnames[MFI_markers]), 1 , paste, collapse="")))
  CV_metaclusters <- matrix(0,
                             nrow = length(files),
                             ncol = max(as.numeric(fsom$metaclustering)) * length(MFI_markers),
                             dimnames = list(names(files),
                                             apply(expand.grid(paste0("MC",1:max(as.numeric(fsom$metaclustering)),"_CV_"), 
                                                               fsom$FlowSOM$prettyColnames[MFI_markers]), 1 , paste, collapse="")))
  
  # abbundance inlezen in matrix (met volledig aantal cellen)-------------------
  for(ft in names(files)){
    if (verbose) { message("Mapping ", ft) }
    ff <- read.FCS(files[ft])
    #NB ingevoegd voor Londen
    #ff@parameters@data$desc[6] <- "CD34 Cy5-5"
    #ff@parameters@data$desc[9] <- "HLA-DR"
    #tot hier
    fsom_tmp <- NewData(fsom$FlowSOM, ff)
    if (verbose) { message("New data ", ft) }
    
    t <- table(fsom_tmp$map$mapping[,1])
    counts[ft, paste0("C",names(t))] <- t
    
    MFI_clusters[ft, ] <- as.vector(GetMFIs(fsom_tmp)[, MFI_markers])
    MFI_metaclusters[ft, ] <- as.vector(MetaclusterMFIs(list(FlowSOM = fsom_tmp,
                                                             metaclustering = fsom$metaclustering))[, MFI_markers])
    
   
    CV_clusters[ft, ] <-  as.vector(GetCVs(fsom_tmp)[, MFI_markers])
    CV_metaclusters[ft, ] <- as.vector(MetaclusterCVs(list(FlowSOM = fsom_tmp,
                                                               metaclustering = fsom$metaclustering))[, MFI_markers])
    
    if (verbose) { message("Mapped ", ft) }
  }
  
  
  # Fill in MFI and CV values for empty clusters with training data values------
  for(ft in names(files)){
   
    na_ids <- which(is.na(MFI_clusters[ft, ]))
    if (verbose) { message("Filling ",length(na_ids)," NA MFI ", ft) }
    for(na_id in na_ids){
      feature <- colnames(MFI_clusters)[na_id]
      cluster <- str_match(feature, "pC([0-9]*)_MFI_(.*)")[,2]
      marker <- str_match(feature, "pC([0-9]*)_MFI_(.*)")[,3]
      channel <- getChannelMarker(ff, marker)[,"name"]
      MFI_clusters[ft, na_id] <- GetMFIs(fsom$FlowSOM)[as.numeric(cluster), 
                                                       channel]
    }
    
    na_ids <- which(is.na(CV_clusters[ft, ]))
    if (verbose) { message("Filling ",length(na_ids)," NA CV ", ft) }
    for(na_id in na_ids){
      feature <- colnames(CV_clusters)[na_id]
      cluster <- str_match(feature, "pC([0-9]*)_CV_(.*)")[,2]
      marker <- str_match(feature, "pC([0-9]*)_CV_(.*)")[,3]
      channel <- getChannelMarker(ff, marker)[,"name"]
      CV_clusters[ft, na_id] <- GetCVs(fsom$FlowSOM)[as.numeric(cluster), 
                                                     channel]
    }
    
    na_ids <- which(is.na(MFI_metaclusters[ft, ]))
    if (verbose) { message("Filling ",length(na_ids)," NA MFI meta ", ft) }
    metaclustermfis <- MetaclusterMFIs(fsom)
    for(na_id in na_ids){
      feature <- colnames(MFI_metaclusters)[na_id]
      cluster <- str_match(feature, "MC([0-9]*)_MFI_(.*)")[,2]
      marker <- str_match(feature, "MC([0-9]*)_MFI_(.*)")[,3]
      channel <- getChannelMarker(ff, marker)[,"name"]
      MFI_metaclusters[ft, na_id] <- metaclustermfis[as.numeric(cluster), 
                                                               channel]
      #CV_metaclusters[ft, na_id] <- MetaclusterCVs(fsom)[as.numeric(cluster), 
      #                                                         channel]
    }
    
    na_ids <- which(is.na(CV_metaclusters[ft, ]))
    if (verbose) { message("Filling ",length(na_ids)," NA CV meta ", ft) }
    metaclustercvs <- MetaclusterCVs(fsom)
    for(na_id in na_ids){
      feature <- colnames(CV_metaclusters)[na_id]
      cluster <- str_match(feature, "MC([0-9]*)_CV_(.*)")[,2]
      marker <- str_match(feature, "MC([0-9]*)_CV_(.*)")[,3]
      channel <- getChannelMarker(ff, marker)[,"name"]

      CV_metaclusters[ft, na_id] <- metaclustercvs[as.numeric(cluster), 
                                                                 channel]
    }
  }

  
  if (verbose) { message("Mapped individual files.")}
  
  pctgs <- t(apply(counts, 1, function(x){x/sum(x)}))
  colnames(pctgs) <- paste0("p", colnames(pctgs))
  meta_pctgs <- t(apply(pctgs, 1, function(x) tapply(x, fsom$metaclustering, sum)))
  colnames(meta_pctgs) <- paste0("MC",colnames(meta_pctgs))
  rel_pctgs <- pctgs
  for(i in seq_len(ncol(pctgs))){
    rel_pctgs[,i] <- rel_pctgs[,i] / meta_pctgs[,fsom$metaclustering[i]] 
  }
  colnames(rel_pctgs) <- paste0("r",colnames(rel_pctgs))
  
  # Replace NAN values with 1 (occurs when a metacluster has only one cluster,
  # and it is empty in one of the files)
  rel_pctgs[is.nan(rel_pctgs)] <- 1 
  
  features <- data.frame(counts,
                         pctgs,
                         meta_pctgs,
                         rel_pctgs,
                         MFI_clusters,
                         MFI_metaclusters,
                         CV_clusters,
                         CV_metaclusters,
                         check.names = FALSE)  
  
  if(!is.null(metadata)) {
    full_data <- left_join(rownames_to_column(features, "Fenotypering_nr"), 
                           metadata, 
                           by = "Fenotypering_nr")
  } else {
    full_data <- rownames_to_column(features, "Fenotypering_nr")
  }
  
  rownames(full_data) <- full_data$Fenotypering_nr
  
  return(full_data)
}




#' Train a FlowSOM and Classification model to predict MDS ---------------------
#' 
#' @param metadata        We assume a column Fenotyperings_nr and MDS
#' @param files_train     Named array with full path to training fcs files and
#'                        phenotype numbers as names
#' @param cols_to_use     Column names or ids from the fcs files to use in the
#'                        FlowSOM clustering
#' @param features_to_use Array with which feature types are passed to the RF,
#'                        options can include "C", "pC", "rC", "MC"
#' 
#' @return Returns a list in which the first item is a FlowSOM result and
#'         the second item a RF model

# Compute FlowSOM for train dataset --------------------------------------------
FlowSOM_training <- function(metadata,
                             files_train,
                             cols_to_use,
                             nCellsPerFile = 1000,
                             gridSize = 15,
                             nClus = 20,
                             seed = NULL,
                             verbose = TRUE){
  
  # Aggregate all train files --------------------------------------------------
  
  if (!is.null(seed)) {set.seed(seed)}
  ff_agg <- AggregateFlowFrames(files_train,
                                cTotal = nCellsPerFile*length(files_train))
  
  if(verbose){message("Aggregated files")}
  
  # FlowSOM --------------------------------------------------------------------
  
  prettyMarkerNames <- ff_agg@parameters@data[,"desc"]
  prettyMarkerNames[is.na(prettyMarkerNames)] <- ff_agg@parameters@data[,"name"][is.na(prettyMarkerNames)]
  names(prettyMarkerNames) <- colnames(ff_agg)
  
  fsom <- FlowSOM(ff_agg,
                  colsToUse = cols_to_use,
                  scale = TRUE,
                  xdim = gridSize, ydim = gridSize,
                  nClus = nClus,
                  seed = seed)
  fsom$FlowSOM$prettyColnames <- prettyMarkerNames
  
  if(verbose){message("Computed FlowSOM")}
  
  # Mapping all files to FlowSOM result ----------------------------------------
  
  full_data <- new_files_to_FlowSOM(fsom,
                                    files_train,
                                    metadata,
                                    verbose = TRUE)
  
  
  
  if(verbose){message("Extracted FlowSOM features")}
  
  return(list(FlowSOM = fsom,
              full_data = full_data))
  
}
# Feature selection ------------------------------------------------------------

FS_byName <- function(data, feature_types){
  cols_for_classification <- NULL
  for (feature_type in feature_types) {
    cols_for_classification <- c(cols_for_classification,
                                 grep(feature_type,
                                      colnames(data)))
  }
  selected_features <- colnames(data)[cols_for_classification]
  return(selected_features)
}

FS_randomforest <- function(forest,
         nFeatures){
  selected_features <- forest$variable.importance %>% 
    sort(decreasing = T) %>% 
    head(nFeatures) %>% 
    names()
  
  return(selected_features)
}

# train classifier--------------------------------------------------------------

RF_training <- function(data,
         features_to_use,
         num.trees = 50000,
         seed = NULL){  
  
  if (!is.null(seed)) {set.seed(seed)}
  
  forest <- ranger::ranger(
    data = data[, c(features_to_use, "MDS_class")],
    dependent.variable.name = "MDS_class",
    num.trees = num.trees,
    importance = "impurity",
    probability = TRUE)
  
  return(forest)
}


SVM_training <- function(data,
         features_to_use,
         seed = NULL){  
  
  if (!is.null(seed)) {set.seed(seed)}
  
  svmmodel <- e1071::svm(data[, features_to_use],
                         data[,"MDS_class"],
                         probability = TRUE)

  return(svmmodel)
}


GMnet_training <- function(data,
                           features_to_use,
                           seed = NULL){  
  
  if (!is.null(seed)) {set.seed(seed)}
  gmnet=glmnet::glmnet(as.matrix(data[, features_to_use]),
                       data[,"MDS_class"],
                       family="binomial")

  return(gmnet)
}




# test classifier --------------------------------------------------------------
model_testing <- function(model,
                          data,
                          classification_method,
                          plot = TRUE){

  if (classification_method == "randomForest") {
    probabilities <- predict(model$classifier,
                     data[ , model$features_to_use])$predictions   
    
    
    #rf_pred <- predict(model$classifier,
    #                         data[ , model$features_to_use])
    #probabilities <- rf_pred$predictions
  } else if (classification_method == "SVM") { 
    svm_pred <- predict(model$classifier,
                        data[ , model$features_to_use],
                        probability = TRUE)
    
    probabilities <- attr(svm_pred, "probabilities")
  } else if (classification_method == "GMnet"){
    gm <- predict(model$classifier,
                  newx = as.matrix(data[ , model$features_to_use]),
                  type = "response",
                  s = min(model$classifier$lambda))
    probabilities <- cbind(1 - gm, gm)
  }
  
  #evaluation of performance sensitivity/specificity/AUC------------------------
  if ("MDS_class" %in% colnames(data)) {
    auc_evaluation <- prediction(probabilities[,2],
                                 data[, "MDS_class"])
    
    cutoff_id <- which.max(performance(auc_evaluation, "spec")@y.values[[1]] + 
                             performance(auc_evaluation, "sens")@y.values[[1]])
    cutoff <- auc_evaluation@cutoffs[[1]][cutoff_id]
    
    if(plot){
      performance_plot <- performance(auc_evaluation, 
                                      measure = "tpr", 
                                      x.measure = "fpr")
      plot(performance_plot, 
           xlab = "population", ylab = "target population", 
           col = "blue", lwd = 2)
      lines(c(0,1),c(0,1),col = "grey")
      
    }
  } else { 
    auc_evaluation = NULL
    cutoff = NULL
  }
  
  return(list(data = data,
              prediction = probabilities,
              evaluation = auc_evaluation,
              best_cutoff = cutoff))
}

# include metadata for files ---------------------------------------------------
read_metadata <- function(file){
  return( readxl::read_xlsx(file) %>% 
            dplyr::filter(Fenotypering_nr > 4000) %>% 
            mutate(DOB = as.Date(as.numeric(DOB), origin="1899-12-30"),
                   Afnamedatum = as.Date(as.numeric(Afnamedatum), origin="1899-12-30"),
                   DOLV = as.Date(as.numeric(DOLV), origin="1899-12-30"),
                   Date_death = as.Date(as.numeric(Date_death), origin="1899-12-30"),
                   MDS = as.factor(MDS),
                   WHO = as.factor(WHO),
                   #`IPSS` = as.factor(`IPSS`),
                   Fenotypering_nr = gsub("ft", "", Fenotypering_nr)) %>% 
            dplyr::filter(Fenotypering_nr > 4000)) # Remove last line !!!
}

get_files_in_dir <- function(directory){
  files_in_dir <- list.files(directory, #pattern=".*fcs$", ###WEGHALEN
                             full.names = TRUE)
  names(files_in_dir) <- gsub(".*([0-9][0-9][0-9][0-9]).*", 
                              "\\1", 
                              files_in_dir)
  return(files_in_dir)
}

# singlet gate------------------------------------------------------------------
is_single <- function(ff, plot = FALSE, ...) {
  fsc_a <- flowCore::exprs(ff)[,"FSC-A"]
  fsc_h <- flowCore::exprs(ff)[,"FSC-H"]
  ratios <- fsc_h / fsc_a
  
  test_points_x <- tapply(fsc_a, cut(fsc_a, 10), mean)
  test_points_y <- tapply(seq_along(fsc_a), cut(fsc_a, 10), function(x){
    diff_fsc_h <- abs(fsc_h[x] - (fsc_a[x] * mean(ratios[x])))
    tr_fsc_h <- mean(diff_fsc_h) + 2 * sd(diff_fsc_h)
    return (mean(fsc_a[x]) * mean(ratios[x]) - tr_fsc_h)
  })
  spl <- splinefun(test_points_x, test_points_y)
  
  selection <- fsc_h > spl(fsc_a)
  return(selection)
}


