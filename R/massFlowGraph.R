#' massFlowGraph
#'
#' A function that permits the user to plot flow frames using the output from
#' `flowPeakDetection()`. The user is prompted to choose which set of frames to
#' plot (all, flagged, or doublet flagged). The chosen flow frames are plotted
#' using the grid functionality in `flowLineGraph()`. The user must ensure the
#' flow directory and file path correspond to the same set of samples for plotting.
#'
#' @param xVariable The fluorescence channel on the x-axis
#' @param flowDir The directory of the gated .fcs data
#' @param filePath The directory of the output .csv file from
#' `flowPeakDetection()` that corresponds to the gated data
#' @param fileName To label the PDF file produced internally by `flowLineGraph()`
#' 
#' @import ggplot2
#' @import tidyverse
#' @import gridExtra
#' @import dplyr
#' @import flowTime
#' 
#' @return description
#' @export
#' 
#' @examples
#' if (interactive()) {
#'   library(BiocFileCache)
#'   original_wd <- getwd()
#'   cachePath <- tempfile(pattern = "bioc_cache_")
#'   dir.create(cachePath, showWarnings = FALSE)
#'   bfc <- BiocFileCache(cachePath, ask = FALSE)
#'   on.exit({
#'     if (normalizePath(cachePath) == normalizePath(getwd())) {
#'       setwd(tempdir())  
#'     }
#'   unlink(cachePath, recursive = TRUE, force = TRUE)
#'   }, add = TRUE)  
#'   gatedFolder <- file.path(cachePath, "fcs_data")
#'   dir.create(gatedFolder, showWarnings = FALSE)
#'   csvFolder <- file.path(cachePath, "csv_data")
#'   dir.create(csvFolder, showWarnings = FALSE)
#'   urls <- c(
#'     "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/gated_data/A1_1.fcs",
#'     "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/gated_data/A1_2.fcs",
#'     "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/gated_data/A1_3.fcs",
#'     "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/gated_data/A1_4.fcs",
#'     "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/gated_data/A1_5.fcs",
#'     "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/gated_data/A1_6.fcs"
#'   )
#'   downloadedFiles <- character(length(urls))
#'   for (i in seq_along(urls)) {
#'     rid <- bfcquery(bfc, urls[i])$rid
#'     if (length(rid) == 0){
#'       rid <- bfcadd(bfc, paste0("file_", i), urls[i])
#'     }
#'     fpath <- bfcrpath(bfc, urls[i])
#'     newFPath <- file.path(gatedFolder, basename(urls[i]))  
#'     file.copy(fpath, newFPath, overwrite = TRUE)
#'     downloadedFiles[i] <- newFPath
#'   }
#'   csv <- "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/analysis/ploidyPeaksOutput.csv"
#'   rid <- bfcquery(bfc, csv, field = "rname")$rid
#'   if (length(rid) == 0) {
#'     rid <- bfcadd(bfc, csv, csv)
#'   }
#'   cachedFile <- bfcrpath(bfc, rid) 
#'   analysisFile <- file.path(csvFolder, basename(csv))
#'   if (!file.exists(analysisFile)) file.copy(cachedFile, analysisFile)
#'   massFlowGraph(
#'     xVariable = "FITC-A",
#'     flowDir = gatedFolder,
#'     filePath = analysisFile,
#'     fileName = "gatedDataRSE"
#'   )
#'}

massFlowGraph <- function(
    xVariable = "FITC-A",
    flowDir = NA,
    filePath = NA,
    fileName = NA
){
  # check for valid file path
  if(!file.exists(filePath)){
    stop("Invalid filepath")
  }
  
  # read in the indicated .csv file 
  data <- read.csv(filePath, header = TRUE)

  # round the final RSE values to 1 decimal
  data$finalRSE <- round(data$finalRSE, 1)
  
  # prompt the user to select what samples they want to plot
  input <- as.numeric(menu(choices = c("Plot All Samples", "Plot Flagged Samples", "Plot Doublet Flagged Samples"), 
                           title = "What would you like to plot?"))
  
  # check if user specifies file name 
  if(!is.na(fileName)){
    output_fileName <- fileName
  }
  
  if(input == 1){
    # check if user specifies file name, give generic name if they don't
    if(is.na(fileName)){
      output_fileName <- "allSamples"
      
    }
  }else if(input == 2){
    # filter data to flagged samples
    data <- data[data$investigate == 1, ]
    
    # check if user specifies file name, give generic name if they don't
    if(is.na(fileName)){
      output_fileName <- "flaggedSamples"
    }
  }else if(input == 3){
    # filter data to doublet flagged samples
    data <- data[data$doublet == 1, ]
    
    # check if user specifies file name, give generic name if they don't
    if(is.na(fileName)){
      output_fileName <- "doubletSamples"
    }
  } else if(input == 0){
    message <- "Exiting the mass plot function"
  }
  
  # get final RSE values for selected samples
  annote <- append(data$finalRSE, "RSE", after = 0)
  
  # if G1_2 peaks are present
  if("G1_2" %in% colnames(data) & (input == 1 | input == 2 | input == 3)){
    flowLineGraph(
      flowDir = flowDir,
      flowSamples = data$Sample,
      xVariable = xVariable,
      grid = TRUE,
      vertLine1 = data$G1_1,
      vertLine2 = data$G1_2,
      annotations = annote,
      fileName = output_fileName
    )
    message <- "Plotting of samples have been completed."
  } else if(input == 1 | input == 2 | input == 3){
    flowLineGraph(
      flowDir = flowDir,
      flowSamples = data$Sample,
      xVariable = xVariable,
      grid = TRUE,
      vertLine1 = data$G1_1,
      annotations = annote,
      fileName = output_fileName
    )
    message <- "Plotting of samples have been completed."
  }
  return(message)
}
