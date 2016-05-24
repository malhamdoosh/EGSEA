#Ensemble of Gene Set Enrichment Analyses 
# 
# Author: Monther Alhamdoosh, E:m.hamdoosh@gmail.com

#' The EGSEAResults class
#'
#' This class contains the results of an EGSEA analysis
#'
#' @slot results list, EGSEA analysis results
#' @slot contrasts character, the analysis contrasts 
#' @slot sampleSize numeric, number of samples
#' @slot gs.annots list, the gene set collection annotation index
#' @slot baseMethods character, vector of base GSE methods
#' @slot combineMethod character, the p-value combining method
#' @slot sort.by character, the results ordering argument
#' @slot symbolsMap data.frame, the mapping between Entrez IDs and Gene Symbols
#' @slot logFC matrix, the logFC matrix of contrasts
#' @slot report logical, whether the report was generated
#' @slot report.dir character, the directory of the EGSEA HTML report
#' 
#' 
#' @name EGSEAResults 
#' @rdname EGSEAResults
#' @aliases EGSEAResults-class
#' @exportClass EGSEAResults

EGSEAResults <- setClass(
            "EGSEAResults",
            
            slots = c(results = "list",
                    contrasts = "character",
                    sampleSize = "numeric",
                    gs.annots = "list",
                    baseMethods = "character",
                    combineMethod = "character",
                    sort.by = "character",
                    symbolsMap = "data.frame",
                    logFC = "matrix",
                    report = "logical",
                    report.dir = "character"),
                      
            prototype = list(results = list(),
                    contrasts = "",
                    sampleSize = 0,
                    gs.annots = list(),
                    baseMethods = c(), 
                    combineMethod = "",
                    sort.by = "",
                    symbolsMap = data.frame(),
                    logFC = matrix(),
                    report = TRUE,
                    report.dir = "./")           
        )
            
   
        
setGeneric(name="addEGSEAResult",
            def = function(object, label, egseaResults){
                standardGeneric("addEGSEAResult")
            }
        )
        
setMethod(f = "addEGSEAResult",
                signature="EGSEAResults",
                definition = function(object, label, egseaResults){
                    object@results[[label]] = egseaResults
                    return(object)
                }
                    
          )
          
