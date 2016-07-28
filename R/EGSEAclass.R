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
#' @rdname EGSEAResults-methods
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
                    symbolsMap = "ANY",
                    logFC = "matrix",
                    report = "logical",
                    report.dir = "character"),                      
            prototype = list(results = list(),
                    contrasts = "",
                    sampleSize = 0,
                    gs.annots = list(),
                    baseMethods = c(), 
                    combineMethod = "fisher",
                    sort.by = "p.adj",
                    symbolsMap = data.frame(),
                    logFC = matrix(),
                    report = TRUE,
                    report.dir = "./")           
        )
     
      
#' @title Extract a slot from an object of class EGSEAResults
#' @description Extract a slot from an object of class EGSEAResults
#' 
#' @param x EGSEAResults object, the analysis result object from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}.
#' @param name character, the slot name
#' 
#' @export
#' @return the selected slot
#' 
#' 
#' @aliases $,EGSEAResults-method
#' @rdname EGSEAResults-methods
#' 
#' @examples
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' class(gsa)
#' print(gsa$baseMethods)

setMethod("$", "EGSEAResults",
        function(x, name){           
            slot(x, name)
        }
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
     
setGeneric(name="addSymbolsMap",
          def = function(object, symbolsMap){
              standardGeneric("addSymbolsMap")
          }
  )
  
  setMethod(f = "addSymbolsMap",
          signature="EGSEAResults",
          definition = function(object, symbolsMap){
              object@symbolsMap = symbolsMap
              return(object)
          }
  
  )
  
  
#' @title Table of Top Gene Sets from an EGSEA Analysis
#' @description Extract a table of the top-ranked gene sets from an EGSEA 
#' analysis.
#' 
#' @param object EGSEAResults object, the analysis result object from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' @param contrast contrast column number or column name specifying which 
#' contrast is of interest.
#' if contrast = 0 or "comparison" and the number of contrasts greater than 1, 
#' the comparative gene sets are 
#' retruned. 
#' @param gs.label the number or label of the gene set collection of interest.
#' @param sort.by character, determines how to order the analysis results in 
#' the stats table. 
#' The accepted values depend on the function used to generate the EGSEA 
#' results.
#' @param number integer, maximum number of gene sets to list
#' @param names.only logical, whether to display the EGSEA statistics or not. 
#' @param verbose logical, whether to print out progress messages and warnings. 
#' 
#' @export
#' @return 
#' A dataframe of top gene sets with the calculated statistics for each if 
#' names.only = FALSE.
#' 
#' @name topSets
#' @aliases topSets,EGSEAResults-method
#' @rdname EGSEAResults-methods
#' 
#' @examples
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' class(gsa)
#' topSets(gsa,  gs.label="kegg",contrast=1, number = 10)
#' topSets(gsa,  gs.label=1, contrast=1, sort.by="ora", number = 10, 
#' names.only=FALSE)
#' topSets(gsa,  gs.label="kegg",contrast=0, number = 10)
#' 
  
  setGeneric(name="topSets",
          def = function(object, gs.label=1, contrast=1, sort.by=NULL, 
                  number = 10, names.only=TRUE, verbose = TRUE){
              standardGeneric("topSets")
          }
  )
  
  setMethod(f = "topSets",
          signature(object = "EGSEAResults"),
          definition = function(object, gs.label=1, contrast=1, sort.by=NULL, 
                  number = 10, names.only=TRUE, verbose = TRUE){
              tryCatch({              
                          if (is.numeric(contrast))
                              if (contrast != 0)
                                  contrast = object@contrasts[contrast]
                              else
                                  contrast = "comparison"
                          if (verbose)
                              cat(paste0("Extracting the top gene sets of the collection \n",
                                              object@gs.annots[[gs.label]]$name, " for the contrast ",
                                              contrast, "\n Sorted by ", 
                                              ifelse(is.null(sort.by), object@sort.by, sort.by)
                                              , "\n"))
                          if (tolower(contrast) == "comparison") 
                              top.gs = object@results[[gs.label]][["comparison"]][["test.results"]]
                          else
                              top.gs = object@results[[gs.label]][["test.results"]][[contrast]]
                          if (! is.null(sort.by)){
                              top.gs = top.gs[order(top.gs[,sort.by],
                                              decreasing=(sort.by == "Significance")), ]
                          }   
                          top.gs = cbind(Rank=seq(1, nrow(top.gs)), top.gs)
                          
                          number = ifelse(number <= nrow(top.gs), number, 
                                  nrow(top.gs))
                          #print(names(top.gs))
                          top.gs = top.gs[1:number, ]
                          if (names.only)
                              return(rownames(top.gs))
                          else{
                              top.gs = as.data.frame(top.gs)
                              top.gs[, "Direction"] = 
                                      as.character(lapply(as.numeric(top.gs[, "Direction"]), 
                                                      function(x) if (x > 0) "Up" else if (x < 0) 
                                                              "Down" else "No Change"))
                              return(top.gs)
                          }
                      }, 
                      error = function(e){
                          cat(paste0("ERROR: topSets(...) encountered an error:\n", e ))
                      })
              return(NULL)
          }
  
  )
  
  
#' @title Show the content of the EGSEAResults object
#' @description Show the parameters of the EGSEAResults object
#' 
#' @inheritParams object EGSEAResults object, the analysis result object from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' 
#' @export
#' @return nothing. 
#' 
#' 
#' @aliases show,EGSEAResults-method
#' @rdname EGSEAResults-methods
#'  
#' @examples
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' class(gsa)
#' show(gsa)
  
setMethod(f = "show",
          signature(object = "EGSEAResults"),
          definition = function(object){
              cat("An object of class \"EGSEAResults\"\n")
              if (length(object@gs.annots) == 0){
                  cat("\tNo results are bound to this object.")
                  return()
              }
              cat(paste0("\tTotal number of genes: ", 
                              length(object@gs.annots[[1]]$featureIDs), "\n"))
              cat(paste0("\tTotal number of samples: ", object@sampleSize, "\n"))
              cat(paste0("\tContrasts: ", paste(object@contrasts, collapse=", "), "\n"))
              cat(paste0("\tBase GSE methods: ", 
                              paste(object@baseMethods, collapse=", "), "\n"))
              cat(paste0("\tP-values combining method: ", object@combineMethod, "\n"))
              cat(paste0("\tSorting statistic: ", object@sort.by, "\n"))
              cat(paste0("\tOrganism: ", object@gs.annots[[1]]$species, "\n"))
              cat(paste0("\tHTML report generated: ", ifelse(object@report, "Yes", "No"), "\n"))
              if (object@report)
                  cat(paste0("\tHTML report directory: ", object@report.dir, "\n"))
              cat("\tTested gene set collections: \n")
              
              for (gs.annot in object@gs.annots){
                  cat("\t\t")
                  #cat(class(gs.annot))
                  #summary(gs.annot)
                  cat(paste0(gs.annot@name, " (", gs.annot@label , "): ",
                                length(gs.annot@idx), " gene sets"))
                  cat("\n")
              }
              cat("Use summary(object) and topSets(object, ...) to explore this object.\n")
          }
)
  
  
#' @title Summary of the EGSEAResults object
#' @description Display a brief summary of the analysis of the EGSEAResults object
#' 
#' @inheritParams object EGSEAResults object, the analysis result object 
#' from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' 
#' @export
#' @return  nothing. 
#' 
#' 
#' @aliases summary,EGSEAResults-method
#' @rdname EGSEAResults-methods
#' 
#' @examples
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' class(gsa)
#' summary(gsa)
  
#setGeneric(name="summary",
#        def = function(object){
#            standardGeneric("summary")
#        }
#)
  
  setMethod(f = "summary",
          signature(object = "EGSEAResults"),
          definition = function(object){
              # for each collection and contrast, print top 10 + comparison if exists  
                if (length(object@gs.annots) == 0){
                    cat("\tNo results are bound to this object.")
                    return()
                }
              for (label in names(object@results)){
                  cat(paste0("**** Top 10 gene sets in the ", 
                                  object@gs.annots[[label]]$name, 
                                  " collection **** \n"))
                  for (contrast in 1:length(object@contrasts)){     
                      cat(paste0("** Contrast ", object@contrasts[contrast], " **\n"))
                      t = topSets(object, label, contrast, verbose=FALSE)
                      for (i in 1:length(t)){
                          cat(t[i])
                          if (i %% 2 == 1)
                              cat(" | ")
                          else
                              cat("\n")
                      }
                      cat("\n")
                  }
                  if ("comparison" %in% names(object@results[[label]])){
                      cat(paste0("** Comparison analysis ** \n"))
                      t = topSets(object, label, "comparison", verbose=FALSE)
                      for (i in 1:length(t)){
                        cat(t[i])
                        if (i %% 2 == 1)
                            cat(" | ")
                        else
                            cat("\n")
                      }
                      cat("\n")
                  }
              }
          }
  
  )
  
#' @title Plot a heatmap for a given gene set
#' @description Generate a heatmap of fold changes for a selected gene set
#' 
#' @inheritParams object EGSEAResults object, the analysis result object from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' @param gene.set character, the name of the gene set. 
#' See the output of \code{\link{topSets}}.
#' @inheritParams gs.label character or numeric, the index/name of the gene set collection.
#' See names(object@@gs.annots) for valid values.
#' @inheritParams contrast character or numeric, the index/name of the contrast or 0/"comparison". 
#' See object@@contrasts for valid values.
#' @param file.name character, the prefix of the output file name. 
#' @param format character, takes "pdf" or "png".
#' @inheritParams verbose logical, whether to print out progress messages and warnings.
#' 
#' @export
#' @return Heatmap plot. 
#' 
#' 
#' @name plotHeatmap
#' @aliases plotHeatmap,EGSEAResults-method
#' @rdname EGSEAResults-methods
#' 
#' @examples
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' class(gsa)
#' plotHeatmap(gsa, gs.label="kegg", "Asthma")
#' plotHeatmap(gsa, gs.label="kegg", "Asthma", contrast = "comparison", 
#' file.name = "asthma.hm.cmp")
  
  setGeneric(name="plotHeatmap",
          def = function(object, gene.set, gs.label=1, contrast=1, file.name="heatmap", 
                  format = "pdf", verbose=TRUE){
              standardGeneric("plotHeatmap")
          }
  )
  
  setMethod(f = "plotHeatmap",
          signature="EGSEAResults",
          definition = function(object, gene.set, gs.label=1, contrast=1, 
                  file.name="heatmap", format = "pdf", verbose=TRUE){            
              tryCatch({            
                          if (is.numeric(contrast))
                              if (contrast != 0)
                                  contrast = object@contrasts[contrast]
                              else
                                  contrast = "comparison"
                          if (verbose)
                              cat(paste0("Generating heatmap for ", gene.set, 
                                              " from the collection \n",
                                              object@gs.annots[[gs.label]]$name, " and for the contrast ",
                                              contrast, "\n"))
                          if (contrast %in% object@contrasts){
                              suppressWarnings(generateHeatMap(gene.set, object@gs.annots[[gs.label]], 
                                              object@logFC[, contrast], 
                                              object@symbolsMap, file.name,
                                              format, print.csv = TRUE))
                          }else if (tolower(contrast) == "comparison"){
                              suppressWarnings(generateHeatMap(gene.set, object@gs.annots[[gs.label]], 
                                              object@logFC, 
                                              object@symbolsMap, file.name,
                                              format, print.csv = TRUE))
                          }else{
                              stop("Unrecognized contrast value. 
                                              Use one of the object@contrasts or a numeric value.")
                          }
                      }, 
                      error = function(e){
                          cat(paste0("ERROR: plotHeatmap(...) encountered an error:\n", e ))
                      })            
          }
  )
  
#' @title Plot a pathway map for a given KEGG pathway
#' @description Generate a visual map for a selected KEGG pathway
#' 
#' @inheritParams object EGSEAResults object, the analysis result object from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' @inheritParams gene.set character, the name of the pathway. 
#' See the output of \code{\link{topSets}}.
#' @inheritParams gs.label character or numeric, the index/name of the KEGG pathways collection.
#' See names(object@@gs.annots)[grep("^kegg", names(object@@gs.annots))] for valid values.
#' @inheritParams contrast character or numeric, the index/name of the contrast or 0/"comparison". 
#' See object@@contrasts for valid values.
#' @inheritParams file.name character, the name of the output file without an extension.
#' @inheritParams verbose logical, whether to print out progress messages and warnings.
#' 
#' 
#' @export
#' @return Pathway map plot. 
#' 
#' 
#' @name plotPathway
#' @aliases plotPathway,EGSEAResults-method
#' @rdname EGSEAResults-methods
#' 
#' @examples
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' class(gsa)
#' plotPathway(gsa, gs.label="kegg", "Asthma")
#' plotPathway(gsa, gs.label="kegg", "Asthma", contrast="comparison", 
#' file.name = "asthma.map.cmp")
  
  setGeneric(name="plotPathway",
          def = function(object, gene.set, gs.label=1, contrast=1, file.name="pathway"
                  , verbose=TRUE){
              standardGeneric("plotPathway")
          }
  )
  
  setMethod(f = "plotPathway",
          signature="EGSEAResults",
          definition = function(object, gene.set, gs.label=1, contrast=1, 
                  file.name="pathway", verbose=TRUE){            
              tryCatch({               
                          if (is.numeric(gs.label))
                              gs.label = names(object@gs.annots)[gs.label]
                          stopifnot(length(grep("^kegg", tolower(gs.label))) == 1)
                          if (is.numeric(contrast))
                              if (contrast != 0)
                                  contrast = object@contrasts[contrast]
                              else
                                  contrast = "comparison"
                          if (verbose)
                              cat(paste0("Generating pathway map for ", gene.set, 
                                              " from the collection \n",
                                              object@gs.annots[[gs.label]]$name, " and for the contrast ",
                                              contrast, "\n"))
                          if (contrast %in% object@contrasts){
                              suppressWarnings(generatePathway(gene.set, object@gs.annots[[gs.label]], 
                                              object@logFC[, contrast], 
                                              file.name = file.name))
                          }else if (tolower(contrast) == "comparison"){
                              suppressWarnings(generateComparisonPathway(gene.set, 
                                              object@gs.annots[[gs.label]], 
                                              object@logFC, 
                                              file.name = file.name))
                          }else{
                              stop("Unrecognized contrast value. 
                                              Use one of the object@contrasts or a numeric value.")
                          }
                      }, 
                      error = function(e){
                          cat(paste0("ERROR: plotPathway(...) encountered an error:\n", e ))
                      })            
          }
  )
  
#' @title Plot a multi-dimensional scaling (MDS) plot for the gene set rankings
#' @description Generate a multi-dimensional scaling (MDS) plot for the gene set rankings
#' of the base GSE methods
#' 
#' @inheritParams object EGSEAResults object, the analysis result object from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' @inheritParams gs.label character or numeric, the index/name of the KEGG pathways collection.
#' See names(object@@gs.annots)[grep("^kegg", names(object@@gs.annots))] for valid values.
#' @inheritParams contrast character or numeric, the index/name of the contrast or 0/"comparison". 
#' See object@@contrasts for valid values.
#' @inheritParams file.name character, the name of the output file without an extension.
#' @inheritParams format character, takes "pdf" or "png".
#' @inheritParams verbose logical, whether to print out progress messages and warnings. 
#' 
#' @export
#' @return MDS plot
#' 
#' 
#' @aliases plotMethods,EGSEAResults-method
#' @rdname EGSEAResults-methods
#'  
#' @examples
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' class(gsa)
#' plotMethods(gsa)
  
setGeneric(name="plotMethods",
          def = function(object, gs.label=1, contrast=1, 
                  file.name="methods.mds", format = "pdf",
                  verbose = TRUE){
              standardGeneric("plotMethods")
          }
  )
  
  setMethod(f = "plotMethods",
          signature="EGSEAResults",
          definition = function(object, gs.label=1, contrast=1, 
                  file.name="methods.mds", format = "pdf",
                  verbose = TRUE){                
              tryCatch({         
                          if (is.numeric(contrast))
                              if (contrast != 0)
                                  contrast = object@contrasts[contrast]
                              else
                                  contrast = "comparison"
                          if (tolower(contrast) == "comparison")
                              stop("plotMethods(...) is not supported for the comparison analysis")
                          if (length(object@baseMethods) < 2){
                              stop("plotMethods(...) requires at least two base methods.")
                          }
                          if (verbose)
                              cat(paste0("Generating MDS plot for the collection \n",
                                              object@gs.annots[[gs.label]]$name, " and for the contrast ",
                                              contrast, "\n"))
                          capture.output(generateMDSMethodsPlot(
                                          object@results[[gs.label]][["test.results"]][[contrast]], 
                                          object@baseMethods, 
                                          file.name, format))
                      }, 
                      error = function(e){
                          cat(paste0("ERROR: plotMethods(...) encountered an error:\n", e ))
                      }
              )
              
          }
  )
  
  
#' @title Generate a Summary plot for EGSEA analysis
#' @description Generate a Summary plot for EGSEA analysis
#' 
#' @inheritParams object EGSEAResults object, the analysis result object from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' @inheritParams gs.label character or numeric, the index/name of the KEGG pathways collection.
#' See names(object@@gs.annots)[grep("^kegg", names(object@@gs.annots))] for valid values.
#' @inheritParams contrast character or numeric, the index/name of the contrast. To compare two
#' contrasts, pass their indexes/names.
#' See object@@contrasts for valid values.
#' @inheritParams file.name character, the name of the heatmap file without an extension.
#' @inheritParams format character, takes "pdf" or "png".
#' @param x.axis character, the x-axis of the summary plot. All the 
#' values accepted by the 
#' \strong{sort.by} parameter can be used. Default x.axis="p.value".
#' @param x.cutoff numeric, cut-off threshold to filter the gene sets of 
#' the summary plots
#' based on the values of the \strong{x.axis}. Default 
#' x.cutoff=NULL.
#' @inheritParams verbose logical, whether to print out progress messages and warnings. 
#' 
#' @export
#' @return Summary plot 
#' 
#' 
#' @name plotSummary
#' @aliases plotSummary,EGSEAResults-method
#' @rdname EGSEAResults-methods
#' 
#' @examples
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' class(gsa)
#' plotSummary(gsa)
#' plotSummary(gsa, contrast=c(1,2), file.name = "summary.cmp")
#' 
  setGeneric(name="plotSummary",
          def = function(object, gs.label=1, contrast=1, 
                  file.name="summary", format = "pdf",
                  x.axis = "p.adj", x.cutoff = NULL,
                  verbose = TRUE){
              standardGeneric("plotSummary")
          }
  )
  
  setMethod(f = "plotSummary",
          signature="EGSEAResults",
          definition = function(object, gs.label=1, contrast=1, 
                  file.name="summary", format = "pdf",
                  x.axis = "p.adj", x.cutoff = NULL,
                  verbose = TRUE){                
              tryCatch({         
                          if (is.numeric(contrast))                            
                              contrast = object@contrasts[contrast]                                  
                          if (is.null(x.cutoff)){
                              if (x.axis %in% c("p.value", "p.adj"))
                                  x.cutoff = 1     
                              else
                                  x.cutoff = 10000
                          }                        
                          if (length(contrast) == 2){
                              if (verbose)
                                  cat(paste0("Generating Summary plots for the collection \n",
                                                  object@gs.annots[[gs.label]]$name, " and for the comparison ",
                                                  paste(contrast, collapse=" vs "), "\n"))                        
                              plot.data = generatePlotData.comparison(
                                      object@results[[gs.label]][["test.results"]][contrast], 
                                      object@results[[gs.label]][["comparison"]][["test.results"]], 
                                      object@gs.annots[[gs.label]], 
                                      x.axis, x.cutoff)
                              if (x.axis %in% c("p.value", "p.adj")){
                                  generateSummaryPlots(plot.data, file.name,
                                          paste0("-log10(p-value) for ", contrast[1]),
                                          paste0("-log10(p-value) for ", contrast[2]),
                                          format = format)
                              }else{
                                  generateSummaryPlots(plot.data, file.name,
                                          paste0("-", x.axis, " for ", contrast[1]),
                                          paste0("-", x.axis, " for ", contrast[2]),
                                          format = format)
                              }
                          }else if (length(contrast) == 1){
                              if (verbose)
                                  cat(paste0("Generating Summary plots for the collection \n",
                                                  object@gs.annots[[gs.label]]$name, " and for the contrast ",
                                                  contrast, "\n"))
                              plot.data = generatePlotData(
                                      object@results[[gs.label]]
                                              [["test.results"]][[contrast]], 
                                      object@gs.annots[[gs.label]], 
                                      x.cutoff, x.axis)
                              if (x.axis %in% c("p.value", "p.adj"))
                                  generateSummaryPlots(plot.data, file.name, 
                                          format = format)
                              else
                                  generateSummaryPlots(plot.data, file.name, 
                                          Xlab = paste0("-", x.axis),
                                          format = format)
                          }else{
                              stop("Wrong number of contrasts. Max is 2.")
                          }                 
                      }, 
                      error = function(e){
                          cat(paste0("ERROR: plotSummary(...) encountered an error:\n", e ))
                      }
              )
              
          }
  )
  

#' @title Plot a GO graph for the GO terms collections 
#' @description Generate a graph of the top significant GO terms in a GO term collection.
#' It could be c5 from MSigDB or Gene Ontolog from the GeneSetDB.
#' 
#' @inheritParams object EGSEAResults object, the analysis result object from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' @inheritParams gs.label character or numeric, the index/name of the KEGG pathways collection.
#' See names(object@@gs.annots)[grep("^kegg", names(object@@gs.annots))] for valid values.
#' @inheritParams contrast character or numeric, the index/name of the contrast or 0/"comparison". 
#' See object@@contrasts for valid values.
#' @inheritParams sort.by
#' @param noSig numeric, number of significant GO terms to be displayed. A number larger than 
#' 5 might not work due to the size of the generated graph. 
#' @inheritParams file.name character, the prefix of the output file name without an extension.
#' @inheritParams format
#' @inheritParams verbose logical, whether to print out progress messages and warnings.
 
#' 
#' 
#' @export
#' @return GO graphs of sginificant BP, MF and CC terms.  
#' 
#' 
#' @name plotGOGraph
#' @aliases plotGOGraph,EGSEAResults-method
#' @rdname EGSEAResults-methods
#' 
#' @examples
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' class(gsa)
#' plotGOGraph(gsa, sort.by="avg.rank")
  
  setGeneric(name="plotGOGraph",
          def = function(object, gs.label="c5", contrast=1, sort.by="p.value",
                  noSig = 5, file.name="c5-top-",
                  format="pdf", verbose=TRUE){
              standardGeneric("plotGOGraph")
          }
  )
  
  setMethod(f = "plotGOGraph",
          signature="EGSEAResults",
          definition = function(object, gs.label="c5", contrast=1, sort.by="p.value",
                  noSig = 5, file.name="c5-top-",
                  format="pdf", verbose=TRUE){            
              tryCatch({               
                          if (is.numeric(gs.label))
                              gs.label = names(object@gs.annots)[gs.label]
                          stopifnot("GOID" %in% colnames(object@gs.annots[[gs.label]]@anno))
                          if (is.numeric(contrast))
                              if (contrast != 0)
                                  contrast = object@contrasts[contrast]
                              else
                                  contrast = "comparison"
                          if (verbose)
                              cat(paste0("Generating GO Graphs for the collection",
                                              object@gs.annots[[gs.label]]$name, 
                                              "\n and for the contrast ",
                                              contrast, " based on the ", sort.by, "\n"))
                        
                          if (contrast %in% object@contrasts){
                              results = object@results[[gs.label]][["test.results"]][[contrast]]
                          }else if (tolower(contrast) == "comparison"){
                              results = object@results[[gs.label]][["comparison"]][["test.results"]]
                          }else{
                              stop("Unrecognized contrast value. 
                                              Use one of the object@contrasts or a numeric value.")
                          }
                          suppressWarnings(
                              generateGOGraphs(
                                  egsea.results = results,
                                  gs.annot = object@gs.annots[[gs.label]], 
                                  sort.by = sort.by,                                          
                                  file.name = file.name, 
                                  noSig = noSig, format = format)
                            )
                      }, 
                      error = function(e){
                          cat(paste0("ERROR: plotPathway(...) encountered an error:\n", e ))
                      })            
          }
  )
  
  
  
  
  
#' @title Display the information of a given gene set name
#' @description Print the details of a given gene set name
#' 
#' @inheritParams object EGSEAResults 
#' @inheritParams gs.label
#' @param set.name character, a vector of gene set names as they appear in \code{\link{topSets}}.
#' 
#' @export
#' @return details of a gene set
#' 
#' 
#' @name showSetByName
#' @aliases showSetByName,EGSEAResults-method
#' @rdname EGSEAResults-methods
#' 
#' @examples
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' class(gsa)
#' showSetByName(gsa, "kegg", "Asthma")
  
setGeneric(name="showSetByName",
          def = function(object, gs.label = 1, set.name){
              standardGeneric("showSetByName")
          }
  )
  
setMethod(f = "showSetByName",
          signature(object = "EGSEAResults"),
          definition = function(object, gs.label = 1, set.name){
              if (is.numeric(gs.label))
                  gs.label = names(object@gs.annots)[gs.label]
              stopifnot(set.name %in% names(object@gs.annots[[gs.label]]@original))
              cols = colnames(object@gs.annots[[gs.label]]@anno)
              for (x in set.name){
                  for (col in cols){
                      cat(paste0(col, ": ", object@gs.annots[[gs.label]]@anno[x, col], "\n"))
                  }
                  cat("\n")
              }
          }
  
  )
  
  
  
#' @title Display the information of a given gene set ID
#' @description Print the details of a given gene set ID
#' 
#' @inheritParams object EGSEAResults
#' @inheritParams gs.label
#' @param id character, a vector of gene set IDs as they appears in the 
#' \code{\link{plotSummary}}.
#' 
#' @export
#' @return  details of a gene set
#' 
#' 
#' @name showSetByID
#' @aliases showSetByID,EGSEAResults-method
#' @rdname EGSEAResults-methods
#'  
#' @examples
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' class(gsa)
#' showSetByID(gsa, "kegg", "hsa04060")
  
  setGeneric(name="showSetByID",
          def = function(object, gs.label = 1, id){
              standardGeneric("showSetByID")
          }
  )
  
  setMethod(f = "showSetByID",
          signature(object = "EGSEAResults"),
          definition = function(object, gs.label = 1, id){
              if (is.numeric(gs.label))
                  gs.label = names(object@gs.annots)[gs.label]
              stopifnot("ID" %in% colnames(object@gs.annots[[gs.label]]@anno))
              stopifnot(id %in% object@gs.annots[[gs.label]]@anno$ID)            
              name = which(object@gs.annots[[gs.label]]@anno$ID %in% id)
              cols = colnames(object@gs.annots[[gs.label]]@anno)
              for (x in name){
                  for (col in cols){
                      cat(paste0(col, ": ", object@gs.annots[[gs.label]]@anno[x, col], "\n"))
                  }
                  cat("\n")
              }
          }
  
  )