# Libraries ---------------------------------------------------------------------------------

library('tidyverse')
library('magrittr')
library('flowCore')
library('flowAI')
library('FlowSOM')
library('flowCore')
library('flowDensity')

# Upload files -----------------------------------------------------------------
dir <- ""
files_of_interest <- list.files(dir, pattern= "")

# Quality control

file <- files_of_interest[1]
for (file in files_of_interest) {
  ff_o <- read.FCS(file.path(dir,file))
  ff_qc <- flow_auto_qc(ff_o, fcs_QC = FALSE, fcs_highQ = "_qc")
}


qc_mini <- file.path("qc/QCmini.txt")
pctg <- read.delim(qc_mini, check.names = FALSE) %>% 
  column_to_rownames("Name file") %>% 
  select(contains("%"))

pdf("flowAI_results.pdf", 
    width = 10,
    height = 20)
pheatmap::pheatmap(pctg,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE)
dev.off()

# Doublet/outlier removal, transformation, uniform flowframes and scaling ------
dir_qc <- "resultsQC"
files_of_interest <- list.files(dir_qc , pattern= "")

  # Upload fcs file for creating uniform flowframes ------------------------------
  ff_13 <- read.FCS(file.path(dir_qc, ""))


for(file in files){
  message(file)
  
  # outlier and doublet removal ------------------------------------------------
  ff <- read.FCS(file.path(dir_qc, file))
  colnames(ff) <- gsub(".*:", "", colnames(ff))
  Result<- flowCore::filter(ff, rectangleGate("FSC-A"=c(25000, 262143)))
  ff <- Subset(ff,Result)
  Result <- flowCore::filter(ff, rectangleGate("SSC-A"=c(0, 262143)))
  ff <- Subset(ff,Result)
  singlets <- is_single(ff)
  ff <- ff[singlets]
  
  #create compensation matrix --------------------------------------------------
  if(!is.null(ff@description$INFINISPILL) && 
     ff@description$INFINISPILL != "NA" &&
     nchar(ff@description$INFINISPILL) > 100) {
    comp_tmp <- strsplit(ff@description$INFINISPILL, ",")[[1]]
    comp <- matrix(as.numeric(comp_tmp[11:74]), ncol = 8, byrow = TRUE)
    colnames(comp) <- gsub(":.*", "", comp_tmp[3:10])
    colnames(comp) <- gsub("", "", colnames(comp))
    ff_c <- compensate(ff, comp)
    
    # Transform, parameters, create uniform flowframes -------------------------    
    channels_to_transform <- colnames(comp)
    ff_t <- transform(ff_c, transformList(channels_to_transform, arcsinhTransform(a = 0, b = 1/150, c = 0)))
    ff_t_new <- flowFrame(exprs = flowCore::exprs(ff_t)[,colnames(ff_13)],
                          parameters = ff_13@parameters,
                          description = ff_13@description)
    
    flowCore::exprs(ff_t_new) <-apply(flowCore::exprs(ff_t_new), 2, function(x){
      return((x - quantile(x, 0.001)) / (quantile(x, 0.999) - quantile(x, 0.001)))
    })
    
   
       suppressWarnings(write.FCS(ff_t_new, file.path(preprocessed_dir, gsub("/", "_", file))))
    
  } else if(!is.null(ff@description$SPILL)){
    comp <- ff@description$SPILL
    ff_c <- compensate(ff, comp)
    
    # Transform, parameters, create uniform flowframes and scale data ----------
    channels_to_transform <- colnames(comp)
    ff_t <- transform(ff_c, transformList(channels_to_transform, arcsinhTransform(a = 0, b = 1/150, c = 0)))
    ff_t_new <- flowFrame(exprs = flowCore::exprs(ff_t)[,colnames(ff_13)],
                          parameters = ff_13@parameters,
                          description = ff_13@description)
    
    flowCore::exprs(ff_t_new) <-apply(flowCore::exprs(ff_t_new), 2, function(x){
      return((x - quantile(x, 0.001)) / (quantile(x, 0.999) - quantile(x, 0.001)))
    })
    
    suppressWarnings(write.FCS(ff_t_new, file.path(preprocessed_dir, gsub("/", "_", file))))
  } else{
    warning("No spillover information for ", file, ", currently skipped.")
  }
}


# Read patient class and select patients----------------------------------------

metadata <- readxl::read_xlsx("metadata.xlsx")
metadata$MDS <- as.factor(metadata$MDS)

selection <- metadata %>% 
  dplyr::filter(Preprocess == 1) %>% 
  dplyr::filter(MDS == 1 | MDS == 2) %>% 
  dplyr::pull(Fenotypering_nr)



# Repeat the following for every tube ------------------------------------------
seed <- 1
plot <- TRUE
cols_to_use <- c(1, 3:7, 9:12)
tubes <- 1:6
outputDir <- "results"
dir.create(outputDir)

for(tube in tubes){
  
  message("Processing tube ",tube)
  
  # Combine subset of cells ----------------------------------------------------
  
  preprocessed_dir <- paste0("test tube ", tube)
  files <- list.files(preprocessed_dir, pattern=".*fcs")
  files <- files[str_match(files,"[0-9][0-9][0-9][0-9]") %in% selection]
  metadata <- files[selection]
  
  set.seed(seed)
  ff_agg <- AggregateFlowFrames(file.path(preprocessed_dir,
                                          files),
                                cTotal = 2000*length(files),
                                writeOutput = FALSE,
                                outputFile = paste0(outputDir,"/Tube_00",tube,"_agg",seed,".fcs"))
  message("  Aggregated ",length(files)," files.")
  
  #  Plot marker expression overview -------------------------------------------
  
  prettyMarkerNames <- ff_agg@parameters@data[,"desc"]
  prettyMarkerNames[is.na(prettyMarkerNames)] <- ff_agg@parameters@data[,"name"][is.na(prettyMarkerNames)]
  names(prettyMarkerNames) <- colnames(ff_agg)
  
  
  if(plot){
    markersToPlot <- colnames(ff_agg)[cols_to_use]
    
    png(file.path(outputDir, paste0("AggregateMarkerOverview_",tube,".png")),
        width = 5000,
        height = 2500)
    par(cex.lab = 2.5, mar = c(4.1,5.1,2.1,2.1))
    layout(matrix(1:12, nrow=3, byrow = TRUE))
    for(marker in markersToPlot){
      print(paste0("Plotting ",prettyMarkerNames[marker]," for the aggregated file."))
      plot(exprs(ff_agg[,c("File_scattered",marker)]),
           pch = ".",
           xlab = "Files",
           ylab = prettyMarkerNames[marker])
    }
    dev.off()
  }
  
  # Flowsom 
  fsom <- FlowSOM(ff_agg,
                  colsToUse = cols_to_use,
                  scale = TRUE,
                  xdim = 15, ydim = 15,
                  nClus = 30,
                  seed = seed)
  fsom$FlowSOM$prettyColnames <- prettyMarkerNames
  fsom$FlowSOM <- BuildMST(fsom$FlowSOM, tSNE=TRUE)
  saveRDS(fsom, file = paste0(outputDir, "/FlowSOM_Tube",tube,
                              "_Markers_",paste(cols_to_use, collapse="_")))
  
  message("  Computed FlowSOM.")
  
  if(plot){
    pdf(file.path(outputDir, paste0("FlowSOM_Tube",tube,".pdf")),
        useDingbats = FALSE)
    
    #Flowsom with legend
    PlotStars(UpdateNodeSize(fsom$FlowSOM, maxNodeSize = 8, reset = TRUE))
    PlotStars(UpdateNodeSize(fsom$FlowSOM, maxNodeSize = 8, reset = TRUE),
              view = "grid")
    
    #Flowsom with metaclustering
    PlotStars(UpdateNodeSize(fsom$FlowSOM, maxNodeSize = 8, reset = TRUE),
              backgroundValues = fsom$metaclustering)
    dev.off()
  }
  
  # abbundance per node matrix 
  counts <- matrix(0,
                   nrow = length(files),
                   ncol = fsom$FlowSOM$map$nNodes,
                   dimnames = list(files,
                                   as.character(1:fsom$FlowSOM$map$nNodes)))
  
  
  # determine abundance for all cells in matrix 
  for(i in seq_along(files)){
    file <- files[i]
    message(file)
    ff <- read.FCS(file.path(preprocessed_dir,file))
    fsom_tmp <- NewData(fsom$FlowSOM, ff)
    
    t <- table(fsom_tmp$map$mapping[,1])
    counts[file, names(t)] <- t
  }
  
  message("  Mapped individual files.")
  
  pctgs <- t(apply(counts, 1, function(x){x/sum(x)}))
  rownames(pctgs) <- gsub(".*([0-9][0-9][0-9][0-9]).*", "\\1", rownames(pctgs))
  pctgs_tmp <- rownames_to_column(data.frame(pctgs), "Fenotypering_nr")
  pctgs_tmp$Fenotypering_nr <- as.factor(pctgs_tmp$Fenotypering_nr)
  metadata$Fenotypering_nr <- as.factor(metadata$Fenotypering_nr)
  pctgs_annotated <- dplyr::left_join(pctgs_tmp, 
                                      metadata, 
                                      by = "Fenotypering_nr")
  pctgs_annotated$Fenotypering_nr <- as.numeric(pctgs_annotated$Fenotypering_nr)
  
  meta_pctgs <- t(apply(pctgs, 1, function(x) tapply(x, fsom$metaclustering, sum)))
  
  rel_pctgs <- pctgs
  
  for(i in seq_len(ncol(pctgs))){
    rel_pctgs[,i] <- rel_pctgs[,i] / meta_pctgs[,fsom$metaclustering[i]] 
  }
  
  save(pctgs_annotated, meta_pctgs, rel_pctgs, file = paste0("pctgs_Tube",tube,
                                                             "_Markers_",paste(cols_to_use, collapse="_"),
                                                             ".Rdata"))
  
  if(plot){
    pdf(paste0("Heatmap_Tube",tube,".pdf"))
    # Heatmap en PCA 
    annotation_row <- data.frame(fileRank = seq_along(files))
    rownames(annotation_row) <- gsub(".*([0-9][0-9][0-9][0-9]).*", "\\1", files)
    
    annotation_col <- data.frame(Metacluster = as.numeric(fsom$metaclustering))
    rownames(annotation_col) <- colnames(pctgs)
      
    pheatmap::pheatmap(meta_pctgs,
                       scale = "column",
                       annotation_col = annotation_col,
                       annotation_row = annotation_row,
                       show_rownames = F,
                       show_colnames = F,
                       main = "Percentages",
                       fontsize = 20)
                       
    pca <- prcomp(scale(meta_pctgs))
    toPlot <- data.frame(File = files,
                         pca_1 = pca$x[,1],
                         pca_2 = pca$x[,2],
                         fileRank = seq_along(files))
    
  ggplot(toPlot) +
      geom_point(aes(x = pca_1, y = pca_2, col = fileRank), size = 4) +
      theme_minimal() +
      theme(text = element_text(size = 12)) +
    scale_color_gradient(low="yellow", high="blue") +
    theme(text = element_text(size=30),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.line = element_line(colour = "black", size = 2),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) 
    
    dev.off()
  }
}

