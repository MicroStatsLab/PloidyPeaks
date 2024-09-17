#' rectGateFlowFrame
#'
#' This function will gate your raw .fcs file and save your gated data, a plot
#' of your gated data, and a plot of the original data
#'
#' @param rawDir The directory of the raw .fcs data
#' @param flowName The name of the .fcs file
#' @param xVariable The fluorescence channel on the x axis
#' @param yVariable The fluorescence channel on the y axis
#' @param xMinValue The lower bound x value for the gate
#' @param xMaxValue The upper bound x value for the gate
#' @param yMinValue The lower bound y value for the gate
#' @param yMaxValue The upper bound y value for the gate
#' @param savePlot A side by side graph comparison the raw data and the gated
#'  data - by default TRUE
#'
#' @return A .fcs of the gated data and plots of gated data
#' 
#' @import patchwork
#' @import tcltk
#' @import flowTime
#' 
#' @export
#'
#' @examples
#' library(BiocFileCache)
#' bfc <- BiocFileCache()
#' url <- "https://ploidypeaksvignette.blob.core.windows.net/ploidypeaksvignettedata/raw_data/A1_1.fcs"
#' subDir <- file.path(bfc@cache, "raw_data")
#' dir.create(subDir, showWarnings = FALSE, recursive = TRUE)
#' targetFile <- file.path(subDir, basename(url))
#' if (!file.exists(targetFile)) file.copy(bfcadd(bfc, url), targetFile)
#' rectGateFlowFrame(
#'  rawDir = targetFile,
#'  flowName = "A1_1.fcs",
#'  xVariable = "FITC-A",
#'  yVariable = "SSC-A",
#'  xMinValue = 50,
#'  xMaxValue = 800,
#'  yMinValue = 50,
#'  yMaxValue = 800,
#'  savePlot = TRUE
#')
#'
rectGateFlowFrame = function(
    rawDir = NA,
    flowName,
    xVariable = "FL1-A",
    yVariable = "SSC-A",
    xMinValue = 10000,
    xMaxValue = 900000,
    yMinValue = 10000,
    yMaxValue = 900000,
    savePlot = TRUE
){
    ##Removing NOTE 'no visible binding for global variable'
    rectGate<-NULL
    
    if(is.na(rawDir)){
        getwd()
        rawDir <- tclvalue(tkchooseDirectory())
    }
    
    #Checking if the folder they selected is empty
    if(purrr::is_empty(rawDir)){
        stop("Your directory is empty")
    }
    
    if(!flowName %in% list.files(rawDir)){
        errorMsg<-paste0(
            "The flow frame ",flowName," is not in the folder, check on the
            spelling of flowName and/or make sure you selected the proper
            folder"
        )
        stop(errorMsg)
        rm(errorMsg)
    }
    
    ##Reading in flow data
    xpectr::suppress_mw(
        flowData <- flowCore::read.FCS(
            paste0(rawDir,"/",flowName),
            transformation=FALSE,
            truncate_max_range=TRUE
        ) 
    )
    
    ##Checking to see if the user input the correct X and Y variable
    if(!xVariable %in% flowData@parameters@data$name){
      stop("Your X variable is not in the dataset")
    }
    
    if(!yVariable %in% flowData@parameters@data$name){
      stop("Your Y variable is not in the dataset")
    }
    
    ##Checking the gating parameters are in the dataset
    if(xMaxValue > max(flowData@exprs[,xVariable])){
      warning("Your xMaxValue exceeds the range of the flow frame. You may want
              to consider a new value")
    }
    
    if(xMinValue < min(flowData@exprs[,xVariable])){
      warning("Your xMinValue is less than the range of the flow frame. You may
              want to consider a new value")
    }
    
    if(yMaxValue > max(flowData@exprs[,yVariable])){
      warning("Your yMaxValue exceeds the range of the flow frame. You may want
              to consider a new value")
    }
    
    if(yMinValue < min(flowData@exprs[,yVariable])){
      warning("Your yMinValue is less than the range of the flow frame. You may
              want to consider a new value")
    }
    
    ##Creating the gate
    autoGate <- paste0(
        'rectGate <- flowCore::rectangleGate(
        filterId=\"Fluorescence Region\",\"',
        xVariable,'\" = c(',xMinValue,',', xMaxValue,'),\"',
        yVariable, '\" = c(',yMinValue,',', yMaxValue,')
        )'
    )
    
    eval(parse(text=autoGate))
    
    ##Subsetting the data that is in the gate
    gatedFlowData <- flowCore::Subset(flowData, rectGate)
    gatedFlowData@description[["GUID"]] <- flowName
    ##Finding the % of cells gated out
    numCellsGatedOut <-  round(
        100 - (
        length(gatedFlowData@exprs[, xVariable]) /
        length(flowData@exprs[, xVariable]))*100, 1
    )
    print( paste0(numCellsGatedOut, "% of the cells were gated out") )
    ##Creating directories to save the data
    setwd(rawDir)
    subDir <- "gated_data"
    dir.create(file.path(dirname(rawDir), subDir), showWarnings=FALSE)
    ##Saving the gated data into a folder
    outFile <- file.path(dirname(rawDir), subDir, flowName)
    flowCore::write.FCS(gatedFlowData, outFile)
    ##If TRUE plots will be created for the user and saved in a folder
    if(savePlot == TRUE){
    
        gatePlotDir <- "plotted_data"
        dir.create(file.path(dirname(rawDir), gatePlotDir), showWarnings=FALSE)
        plotOutFile <- file.path(dirname(rawDir), gatePlotDir)
        
        flowData@description[["GUID"]] <- "Raw data"
        rawDataPlot <- ggcyto::autoplot(
            flowData, xVariable, yVariable, bins=64) + 
            ggcyto::geom_gate(rectGate) + 
            ggcyto::ggcyto_par_set(limits="data")
        gatedFlowData@description[["GUID"]] <- "Gated data"
        gatedDataPlot <- xpectr::suppress_mw(
            ggcyto::autoplot(gatedFlowData, xVariable, yVariable,  bins=64)
            # + ggcyto::ggcyto_par_set(limits="instrument")
        )
        combinedPlot <- ggcyto::as.ggplot(rawDataPlot) + 
            ggcyto::as.ggplot(gatedDataPlot)
    
        gatedFlowData@description[["GUID"]] <- flowName
        png(paste0(plotOutFile, "/", flowName, '.png'), width=600, height=400)
        print(combinedPlot)
        dev.off()
        
    }else{
        flowData@description[["GUID"]] <- "Raw data"
        rawDataPlot <- ggcyto::autoplot(
            flowData, xVariable, yVariable, bins=64) + 
            ggcyto::geom_gate(rectGate) + 
            ggcyto::ggcyto_par_set(limits="data")
        gatedFlowData@description[["GUID"]] <- "Gated data"
        gatedDataPlot <- xpectr::suppress_mw(
            ggcyto::autoplot(gatedFlowData, xVariable, yVariable,  bins=64) + 
                ggcyto::ggcyto_par_set(limits="instrument")
        )
        combinedPlot <- ggcyto::as.ggplot(rawDataPlot) + 
            ggcyto::as.ggplot(gatedDataPlot)
        gatedFlowData@description[["GUID"]] <- flowName
    }
    return(combinedPlot)
}
