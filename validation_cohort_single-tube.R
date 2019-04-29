source("helperFunctions_AI_workflow.R")

# Final diagnostic pipeline ----------------------------------------------------

# Directory in which the control sample is located
# This sample is used to choose uniform channel/marker names
control_files_dir <- ""


# Directory in which the raw fcs files are saved
id <- ""

dir <- paste0("", id)

# Tube numbers of interest
tubes <- c("4")
# Files in the given directory which match the tube numbers
files <- list.files(dir)

files <- sapply(tubes,
                function(x) grep(paste0("_00", x, ".fcs$"), 
                                 files, 
                                 value = TRUE))



# Directory to save the preprocessed fcs files in
preprocessed_dir <- paste0(dir, "_preprocessed")
dir.create(preprocessed_dir)
feature_types <- c(#pC= "^pC[0-9]*", 
                   #rpC = "^rpC[0-9]*", 
                   MC = "^MC[0-9]*")

outputDir <- "Final_model"

# 1. Preprocessing -------------------------------------------------------------
tube <- tubes[1]
for(tube in tubes){

  # Read the right control file for this tube
  control_file <- paste0("", tube, "_qc.fcs")
  ff_control <- read.FCS(file.path(control_files_dir, 
                                   control_file))
  
  #Actual file of interest
  file <- files[tube]
  
  # 1.1 FlowAI -----------------------------------------------------------------
  
  ff <- flow_auto_qc(file.path(dir, file), 
                     fcs_QC = FALSE, 
                     fcs_highQ = FALSE)
  colnames(ff) <- gsub(".*:", "", colnames(ff))
  
  
  # 1.2 Pre-gaten, compenseren, transformeren, gelijke flowframes, scalen ------
  
  # Remove boundary values
  filter_fsc <- flowCore::filter(ff, rectangleGate("FSC-A" = c(25000, 262143)))
  ff <- Subset(ff, filter_fsc)
  filter_ssc <- flowCore::filter(ff, rectangleGate("SSC-A" = c(0, 262143)))
  ff <- Subset(ff, filter_ssc)
  
  # Remove doublets
  singlets <- is_single(ff)
  ff <- ff[singlets, ]
  
  # Compensate
  # Compensation matrix might be saved in 2 ways, either in SPILL or INFINISPILL
  if(!is.null(ff@description$SPILL)){
    
    comp <- ff@description$SPILL
  
  } else if (!is.null(ff@description$INFINISPILL) && 
             ff@description$INFINISPILL != "NA" &&
             nchar(ff@description$INFINISPILL) > 100) {
    
    # Get spillover information from a string to a matrix
    comp_tmp <- strsplit(ff@description$INFINISPILL, ",")[[1]]
    comp <- matrix(as.numeric(comp_tmp[11:74]), ncol = 8, byrow = TRUE)
    colnames(comp) <- gsub(":.*", "", comp_tmp[3:10])
    colnames(comp) <- gsub("CD34 Cy5-5", "PerCP-A", colnames(comp))
    
  } else {
    warning("No spillover information for ", file, ", currently skipped.")
    comp <- NULL
  }
  
  if (!is.null(comp)){
    ff <- compensate(ff, comp)
    
    # Transform
    channels_to_transform <- colnames(comp)
    ff <- transform(ff, transformList(channels_to_transform, 
                                      arcsinhTransform(a = 0, b = 1/150, c = 0)))
    
    # New flowframe to save with the right channel and marker names
    ff <- flowFrame(exprs = flowCore::exprs(ff)[,colnames(ff_control)],
                    parameters = ff_control@parameters,
                    description = ff_control@description)
    
    # Rescale data, with 0.01 and 0.99 quantiles resp set to 0 and 1
    flowCore::exprs(ff) <-apply(flowCore::exprs(ff), 2, function(x){
      return((x - quantile(x, 0.001)) / (quantile(x, 0.999) - quantile(x, 0.001)))
    })
    
    # Write result to new fcs file
    suppressWarnings(write.FCS(ff, 
                               file.path(preprocessed_dir, 
                                         gsub("/", "_", file))))
    
    
  }
} 

# 2. Plotting new file on FlowSOM models ----------------------------------------

version <- ""

features_all <- NULL
for(tube in tubes){
  
  file <- get_files_in_dir(preprocessed_dir)
 
  model_file <- file.path(outputDir,
                          paste0("final_fsommodel",version, tube, ".RDS"))
  fsommodel <- readRDS(model_file)
  features_enr <- new_files_to_FlowSOM(fsommodel$FlowSOM,
                                       files = file,
                                       metadata = NULL,
                                       verbose = TRUE)
  
  features_to_use <- FS_byName(fsommodel$full_data, feature_types)
  feature_ids <- which(colnames(features_enr) %in% features_to_use)
  colnames(features_enr)[feature_ids] <-  
    paste0("T",
          tube,
           "_",
           colnames(features_enr)[feature_ids])
  if (is.null(features_all)) {
    features_all <- features_enr
  } else {
    features_all <- full_join(features_all, 
                              features_enr, 
                              by = "Fenotypering_nr")
  }
  
}





# 3. Performing classification on new file -------------------------------------
final_model <- readRDS("")
prediction <-  model_testing(final_model,
                             data = features_all,
                             classification_method = "randomForest",
                             plot = TRUE)

predicted_label <- c("Control", "MDS")[1 + (prediction$prediction[,2] >= 0.5)]
print(prediction$prediction)
print(predicted_label)

write.csv2(prediction$prediction,paste0(dir,"prediction2.csv"))


