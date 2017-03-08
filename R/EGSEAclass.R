#Ensemble of Gene Set Enrichment Analyses 
# 
# Author: Monther Alhamdoosh, E:m.hamdoosh@gmail.com

#' @title The EGSEAResults class
#'
#' @description The \code{EGSEAResults} class stores the results of an EGSEA analysis. 
#' @details The \code{EGSEAResults} class is used by \code{egsea}, \code{egsea.cnt} and
#' \code{egsea.ora} to store the results of an EGSEA analysis. This helps in mining the 
#' analysis results and generating customized tables and plots.   
#'
#' @slot results list, EGSEA analysis results
#' @slot limmaResults MArrayLM, is a limma linear fit model
#' @slot contr.names character, the contrasts defined in the analysis 
#' @slot contrast double, an N x L matrix indicates the contrasts of the 
#' linear model coefficients for which the test is required. N is the number of 
#' columns of the design matrix and L is number of contrasts. Can be also a vector 
#' of integers that specify the columns of the design matrix. 
#' @slot sampleSize numeric, number of samples
#' @slot gs.annots list, the gene set collection annotation index
#' @slot baseMethods character, vector of base GSE methods
#' @slot baseInfo list, additional information on the base methods (e.g., version).
#' @slot combineMethod character, the p-value combining method
#' @slot sort.by character, the results ordering argument
#' @slot symbolsMap data.frame, the mapping between Entrez IDs and Gene Symbols
#' @slot logFC matrix, the logFC matrix of contrasts
#' @slot logFC.calculated character, indicates whether the logFC was calculated using
#' limma DE analysis.
#' @slot sum.plot.axis character, the x-axis of the summary plot
#' @slot sum.plot.cutoff numeric, the cut-off threshold for the summary plot x-axis
#' @slot report logical, whether the report was generated
#' @slot report.dir character, the directory of the EGSEA HTML report
#' @slot egsea.version character, the version of EGSEA package
#' @slot egseaData.version character, the version of EGSEAdata package
#' 
#' @importClassesFrom limma MArrayLM
#' 
#' @name EGSEAResults 
#' @rdname EGSEAResults-methods
#' @aliases EGSEAResults-class
#' @exportClass EGSEAResults

EGSEAResults <- setClass(
            "EGSEAResults",            
            slots = c(results = "list",
                    limmaResults = "MArrayLM",
                    contr.names = "character",
                    contrast = "ANY",
                    sampleSize = "numeric",
                    gs.annots = "list",
                    baseMethods = "character",
                    baseInfo= "list",
                    combineMethod = "character",
                    sort.by = "character",
                    symbolsMap = "ANY",
                    logFC = "matrix",
                    logFC.calculated = "character",
                    sum.plot.axis = "character",
                    sum.plot.cutoff = "ANY",
                    report = "logical",
                    report.dir = "character",
                    egsea.version = "character",
                    egseaData.version = "character"),                      
            prototype = list(results = list(),
                    limmaResults = new("MArrayLM"),
                    contr.names = "",
                    contrast = c(),
                    sampleSize = 0,
                    gs.annots = list(),
                    baseMethods = c(), 
                    baseInfo = list(),
                    combineMethod = "fisher",
                    sort.by = "p.adj",
                    symbolsMap = data.frame(),
                    logFC = matrix(),
                    logFC.calculated = "No",
                    sum.plot.axis = "p.adj",
                    sum.plot.cutoff = NULL,
                    report = TRUE,
                    report.dir = "./",
                    egsea.version = "Unknown",
                    egseaData.version = "Unknown")           
        )
     
      
#' @title Extract a slot from an object of class EGSEAResults
#' @description The opertator \code{$} extracts a slot from an object of class EGSEAResults.
#' 
#' @param x EGSEAResults object, the analysis result object from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}.
#' @param name character, the slot name
#' 
#' @export
#' @return \code{$} returns the selected slot. 
#' 
#' 
#' @aliases $,EGSEAResults-method
#' @rdname EGSEAResults-methods
#' 
#' @examples
#' # Exampple of EGSEAResults
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' class(gsa)
#' print(gsa$baseMethods)
#' 

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
#' @description \code{topSets} extracts a table of the top-ranked gene sets from an EGSEA 
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
#' the stats table. The accepted values depend on the function used to generate the EGSEA 
#' results.
#' @param number integer, maximum number of gene sets to list
#' @param names.only logical, whether to display the EGSEA statistics or not. 
#' @param verbose logical, whether to print out progress messages and warnings. 
#' 
#' @export
#' @return 
#' \code{topSets} returns a dataframe of top gene sets with the calculated statistics for each if 
#' names.only = FALSE.
#' 
#' @name topSets
#' @aliases topSets,EGSEAResults-method
#' @rdname EGSEAResults-methods
#' 
#' @examples
#' # Example of topSets
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
                                  contrast = object@contr.names[contrast]
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
                              top.gs[, "direction"] = 
                                      as.character(lapply(as.numeric(top.gs[, "direction"]), 
                                                      function(x) if (x > 0) "Up" else if (x < 0) 
                                                              "Down" else "Neutral"))
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
#' @description \code{show} displays the parameters of an EGSEAResults object
#' 
#' @inheritParams object EGSEAResults object, the analysis result object from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' 
#' @export
#' @return \code{show} does not return data.  
#' 
#' 
#' @aliases show,EGSEAResults-method
#' @rdname EGSEAResults-methods
#'  
#' @examples
#' # Example of show
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' class(gsa)
#' show(gsa)
#' 
  
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
              cat(paste0("\tContrasts: ", paste(object@contr.names, collapse=", "), "\n"))
              base.names = names(object@baseInfo)
              base.vers = sapply(base.names, function(x) 
                          as.character(object@baseInfo[[x]]$version))
              base.pkgs = sapply(base.names, function(x) object@baseInfo[[x]]$package)
              baseMethods = paste0(base.names, " (", base.pkgs, ":", base.vers, ")")
              cat(paste0("\tBase GSE methods: ", 
                              paste(baseMethods, collapse=", "), "\n"))
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
                                length(gs.annot@idx), " gene sets - Version: ",
                                gs.annot@version, ", Update date: ", gs.annot@date))
                  cat("\n")
              }
              cat(paste0("\tEGSEA version: ", object@egsea.version, "\n"))
              cat(paste0("\tEGSEAdata version: ", object@egseaData.version, "\n"))
              cat("Use summary(object) and topSets(object, ...) to explore this object.\n")
          }
)
  

  
#' @title Summary of the EGSEAResults object
#' @description \code{summary} displays a brief summary of the analysis results 
#' stored in an EGSEAResults object
#' 
#' @inheritParams object EGSEAResults object, the analysis result object 
#' from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' 
#' @export
#' @return \code{summary} does not return data. 
#' 
#' 
#' @aliases summary,EGSEAResults-method
#' @rdname EGSEAResults-methods
#' 
#' @examples
#' # Example of summary
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' class(gsa)
#' summary(gsa)
#' 
  
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
                  for (contrast in 1:length(object@contr.names)){     
                      cat(paste0("** Contrast ", object@contr.names[contrast], " **\n"))
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
  
#' @title Top tables of the limma differential expression analysis
#' @description \code{limmaTopTable} returns a dataframe of the top table of the 
#' limma analysis for a given contrast. 
#' @details \code{limmaTopTable} output can be understood from \code{limma::topTable}. 
#' 
#' @inheritParams object EGSEAResults object, the analysis result object 
#' from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' @inheritParams contrast character or numeric, the index/name of the contrast 
#' 
#' @export
#' @return  \code{limmaTopTable} returns a dataframe. 
#' 
#' 
#' @aliases limmaTopTable,EGSEAResults-method
#' @rdname EGSEAResults-methods
#' 
#' @examples
#' # Example of limmaTopTable
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' class(gsa)
#' colnames(limmaTopTable(gsa))
#' head(limmaTopTable(gsa))
#' 

setGeneric(name="limmaTopTable",
        def = function(object, contrast = 1){
            standardGeneric("limmaTopTable")
        }
)
  
setMethod(f = "limmaTopTable",
          signature(object = "EGSEAResults"),
          definition = function(object, contrast = 1){
              if (length(object@limmaResults) > 0 ){
                  if (is.numeric(contrast))
                      stopifnot(contrast > 0 && contrast <= length(object@contr.names))
                  else
                      stopifnot(contrast %in% object@contr.names)     
                  t = get.toptables(object@limmaResults, object@contrast)[[contrast]]                  
                  rownames(t) = t[,1]
                  return(t[order(t[, "adj.P.Val"]), ])
              }else{
                  cat("Limma analysis results are not available. \n")
                  cat("Try to re-run EGSEA analysis with keep.limma = TRUE. \n")
                  return(NULL)
              }
          }
)

#' @title HTML report of the EGSEA analysis 
#' @description \code{generateReport} creates an HTML report for the EGSEA analysis that
#' enables users to seamlessly browse the test results.  
#' @details EGSEA report is an interactive HTML report that is generated to
#' enable a swift navigation through the results of an EGSEA analysis. The following pages 
#' are generated for each gene set collection and contrast/comparison: \cr
#' 1. Stats Table page shows the detailed statistics of the EGSEA analysis for the 
#' \code{display.top} gene sets. It shows the EGSEA scores, individual rankings and 
#' additional annotation for each gene set. Hyperlinks to the source of each gene set
#' can be seen in this table when they are available. The "Direction" column shows the regulation
#' direction of a gene set which is calculated based on the \code{logFC}, which is
#' either calculated from the limma differential expression analysis or provided by the user. 
#' The method \code{topSets} can be used to generate custom Stats Table. \cr
#' 2. Heatmaps page shows the heatmaps of the gene fold changes for the gene sets that are
#' presented in the Stats Table page. Red indicates up-regulation
#' while blue indicates down-regulation. Only genes that appear in the input expression/count 
#' matrix are visualized in the heat map. Gene names are coloured based on their 
#' statistical significance in the \code{limma} differential expression analysis. 
#' The "Interpret Results" link below each heat map allows the user to download the 
#' original heat map values along with additional statistics from \code{limma} DE analysis (
#' if available) so that they can be used to perform further analysis in R, e.g., customizing 
#' the heat map visualization. Additional heat maps can be generated and customized
#'  using the method \code{plotHeatmap}. \cr
#' 3. Summary Plots page  shows the methods ranking plot along with the summary plots of 
#' EGSEA analysis. The method plot uses multidimensional scaling (MDS) to visualize the 
#' ranking of individual methods on a given gene set collection. The summary plots are 
#' bubble plots that visualize the distribution of gene sets based on the EGSEA
#' Significance Score and another EGSEA score (default, p-value). 
#' Two summary plots are generated: ranking and directional plots. Each gene set is 
#' reprersented with a bubble which is coloured based on the EGSEA ranking (in ranking
#' plots ) or gene set regulation direction (in directional plots) and sized based on the 
#' gene set cardinality (in ranking plots) or EGSEA Significance score (in directional plots).
#' Since the EGSEA "Significance Score" is proportional to the p-value and the 
#' absolute fold changes, it could be useful to highlight gene sets that
#' have high Significance scores. The blue labels on the summary plot indicate 
#' gene sets that do not appear in the top 10 list of gene sets based on the "sort.by" 
#' argument (black labels) yet they appear in the top 5 list of gene sets based on 
#' the EGSEA "Significance Score". If two contrasts are provided, the rank is calculated 
#' based on the "comparison" analysis results and the "Significance Score" is calculated 
#' as the mean. The method \code{plotSummary} can be used to customize the Summary plots by 
#' changing the x-axis score
#' and filtering bubbles based on the values of the x-axis. The method \code{plotMethods} can be
#' used to generate Method plots. \cr
#' 4. Pathways page shows the KEGG pathways for the gene sets that are presented in the
#' Stats Table of a KEGG gene set collection. The gene fold changes are overlaid on the 
#' pathway maps and coloured based on the gene regulation direction: blue for down-regulation
#' and red for up-regulation. The method \code{plotPathway} can be used to generate
#' additional pathway maps. Note that this page only appears if a KEGG gene set collection
#' is used in the EGSEA analysis. \cr
#' 5. Go Graphs page shows the Gene Ontology graphs for top 5 GO terms in each of 
#' three GO categories: Biological Processes (BP), Molecular Functions (MF), 
#' and Cellular Components (CC). Nodes are coloured based on the default \code{sort.by}
#' score where red indicates high significance and yellow indicates low significance. 
#' The method \code{plotGOGraph} can be used to customize GO graphs by 
#' changing the default sorting score and the number of significance nodes that can be
#' visualized. It is recommended that a small number of nodes is selected. Note that
#' this page only appears if a Gene Ontology gene set collection is used, i.e., for
#' the c5 collection from MSigDB or the gsdbgo collection from GeneSetDB. \cr
#' \cr
#' Finally, the "Interpret Results" hyperlink in the EGSEA report allows the user to download
#' the fold changes and limma analysis results and thus improve the interpretation of the results.
#' 
#' @inheritParams object EGSEAResults object, the analysis result object 
#' from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' @inheritParams number
#' @inheritParams sort.by
#' @param egsea.dir character, directory into which the analysis results are 
#' written out. 
#' @param kegg.dir character, the directory of KEGG pathway data file (.xml) 
#' and image file (.png). 
#' Default kegg.dir=paste0(egsea.dir, "/kegg-dir/").
#' @inheritParams x.axis
#' @inheritParams x.cutoff
#' @param num.threads numeric, number of CPU cores to be used. Default 
#' num.threads=4.
#' @param print.base logical, whether to write out the analysis results of the base methods.
#' Default is False.
#' @inheritParams verbose
#' 
#' @export
#' @return  \code{generateReport} does not return data but creates an HTML report.
#' 
#' 
#' @aliases generateReport,EGSEAResults-method
#' @rdname EGSEAResults-methods
#' 
#' @examples
#' # Example of generateReport
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' # generateReport(gsa)
#'

setGeneric(name="generateReport",
        def = function(object, number = 20, sort.by = NULL, 
                egsea.dir = NULL, kegg.dir = NULL, 
                x.axis = NULL, x.cutoff = NULL, 
                num.threads = 4,
                print.base = FALSE,
                verbose = FALSE){
            standardGeneric("generateReport")
        }
)

setMethod(f = "generateReport",
        signature(object = "EGSEAResults"),
        definition = function(object, number = 20, sort.by = NULL, 
                egsea.dir = NULL, kegg.dir = NULL, 
                x.axis = NULL, x.cutoff = NULL, 
                num.threads = 4,
                print.base = FALSE,
                verbose = FALSE){
            generate.EGSEA.Report(object,  
                    display.top = number, sort.by = sort.by, 
                    egsea.dir = egsea.dir, kegg.dir = kegg.dir, 
                    sum.plot.axis = x.axis, sum.plot.cutoff = x.cutoff, 
                    num.threads = num.threads,
                    print.base = print.base,
                    verbose = verbose)
        }
)


#' @title Results of the limma differential expression analysis
#' @description \code{getlimmaResults} returns the linear model fit produced by
#'  \code{limma::eBayes}. 
#' @details \code{getlimmaResults}'s output can be manipulated using
#'  \code{limma::topTable} and \code{limma::topTreat}. 
#' 
#' @inheritParams object EGSEAResults object, the analysis result object 
#' from  \code{\link{egsea}} or
#' \code{\link{egsea.cnt}}. 
#' 
#' @export
#' @return  \code{getlimmaResults} returns an MArrayLM object.
#' 
#' 
#' @aliases getlimmaResults,EGSEAResults-method
#' @rdname EGSEAResults-methods
#' 
#' @examples
#' # Example of getlimmaResults
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' class(gsa)
#' fit = getlimmaResults(gsa)
#' class(fit)
#' names(fit)
#' 

setGeneric(name="getlimmaResults",
        def = function(object){
            standardGeneric("getlimmaResults")
        }
)

setMethod(f = "getlimmaResults",
        signature(object = "EGSEAResults"),
        definition = function(object){
            if (length(object@limmaResults) > 0 ){               
                return(object@limmaResults)
            }else{
                cat("Limma analysis results are not available. \n")
                cat("Try to re-run EGSEA analysis with keep.limma = TRUE. \n")
                return(NULL)
            }
        }
)
  
#' @title Plot a heatmap for a given gene set
#' @description \code{plotHeatmap} generates a heatmap of fold changes for a selected gene set.
#' @details  \code{plotHeatmap} fold changes are colored based on the \code{fc.colors} and 
#' only genes that appear in the EGSEA analysis are visualized in the heatmap. Gene names 
#' are coloured based on the statistical significance level from limma DE analysis. 
#' 
#' @inheritParams object EGSEAResults object, the analysis result object from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' @param gene.set character, the name of the gene set. 
#' See the output of \code{\link{topSets}}.
#' @inheritParams gs.label character or numeric, the index/name of the gene set collection.
#' See \code{names(object@@gs.annots)} for valid values.
#' @inheritParams contrast character or numeric, the index/name of the contrast or 0/"comparison". 
#' See \code{object@@contr.names} for valid values.
#' @param file.name character, the prefix of the output file name. 
#' @param format character, takes "pdf" or "png".
#' @inheritParams verbose logical, whether to print out progress messages and warnings.
#' @param fc.colors vector, determines the fold change colors of the heatmap. 
#' Three colors of the negative, zero and positive log fold changes,
#' respectively, should be assigned. Default is c( "#67A9CF", "#F7F7F7", "#EF8A62"). These 
#' colors were generated using \code{rev(RColorBrewer::brewer.pal(3, "RdBu"))}
#' @export
#' @return \code{plotHeatmap} does not return data but creates image and CSV files. 
#' 
#' 
#' @name plotHeatmap
#' @aliases plotHeatmap,EGSEAResults-method
#' @rdname EGSEAResults-methods
#' 
#' @examples
#' # Example of plotHeatmap
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' class(gsa)
#' plotHeatmap(gsa, "Asthma", gs.label="kegg")
#' plotHeatmap(gsa, "Asthma", gs.label="kegg", contrast = "comparison", 
#' file.name = "asthma.hm.cmp")
#' 
  
  setGeneric(name="plotHeatmap",
          def = function(object, gene.set, gs.label=1, contrast=1, file.name="heatmap", 
                  format = "pdf",
                  fc.colors= c( "#67A9CF", "#F7F7F7", "#EF8A62"),
                  verbose=TRUE){
              standardGeneric("plotHeatmap")
          }
  )
  
  setMethod(f = "plotHeatmap",
          signature="EGSEAResults",
          definition = function(object, gene.set, gs.label=1, contrast=1, 
                  file.name="heatmap", format = "pdf", 
                  fc.colors= c( "#67A9CF", "#F7F7F7", "#EF8A62"),
                  verbose=TRUE){            
              tryCatch({            
                          if (is.numeric(contrast))
                              if (contrast != 0)
                                  contrast = object@contr.names[contrast]
                              else
                                  contrast = "comparison"
                          if (verbose)
                              cat(paste0("Generating heatmap for ", gene.set, 
                                              " from the collection \n",
                                              object@gs.annots[[gs.label]]$name, " and for the contrast ",
                                              contrast, "\n"))
                          if (contrast %in% object@contr.names){
                              if (length(object@limmaResults) > 0){
                                  t = topTable(object@limmaResults, coef=contrast, 
                                      number=Inf, sort.by="none")
                                  rownames(t) = rownames(object@limmaResults)
                                  limma.tops = list(contrast = t)
                              }else{
                                  cat("WARNING: limma analysis results are not available.")
                                  limma.tops = list()
                              }
                              suppressWarnings(generateHeatMap(gene.set, object@gs.annots[[gs.label]], 
                                              object@logFC[, contrast],
                                              limma.tops,
                                              object@symbolsMap, file.name,
                                              format, print.csv = TRUE,
                                              fc.colors))
                          }else if (tolower(contrast) == "comparison"){
                              if (length(object@limmaResults) > 0){
                                  limma.tops = list()
                                  for (c in object@contr.names){
                                     t = topTable(object@limmaResults, coef=c, 
                                          number=Inf, sort.by="none")
                                     rownames(t) = rownames(object@limmaResults)
                                     limma.tops[[c]] = t
                                  }
                              }else{
                                  cat("WARNING: limma analysis results are not available.")
                                  limma.tops = list()
                              }
                              suppressWarnings(generateHeatMap(gene.set, object@gs.annots[[gs.label]], 
                                              object@logFC, 
                                              limma.tops,
                                              object@symbolsMap, file.name,
                                              format, print.csv = TRUE,
                                              fc.colors))
                          }else{
                              stop("Unrecognized contrast value. 
                                              Use one of the object@contr.names or a numeric value.")
                          }
                      }, 
                      error = function(e){
                          cat(paste0("ERROR: plotHeatmap(...) encountered an error:\n", e ))
                      })            
          }
  )
  

#' @title Plot a summary heatmap for EGSEA analysis 
#' @description \code{plotSummaryHeatmap} generates a summary heatmap for the top n gene 
#' sets of the comparative analysis across multiple contrasts. 
#' @details  \code{plotSummaryHeatmap} creates a summary heatmap for the rankings
#' of top \code{number} gene sets of the comparative analysis across all the contrasts. The
#' \code{show.vals} score can be displayed on the heatmap for each gene set. This can
#' help to identify gene sets that are highly ranked/sgnificant across multiple
#' contrasts. 
#' 
#' @inheritParams object EGSEAResults object, the analysis result object from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' @inheritParams gs.label character or numeric, the index/name of the gene set collection.
#' See names(object@@gs.annots) for valid values.
#' @inheritParams number integer, maximum number of gene sets to list
#' @inheritParams sort.by character
#' @param hm.vals character, determines which EGSEA score values are used to draw the map.
#' Default is NULL which implies using the sort.by score. 
#' @param show.vals character, determines which EGSEA score values are displayed on the map.
#' Default is NULL which does not show anything.
#' @inheritParams file.name character, the prefix of the output file name. 
#' @inheritParams format character, takes "pdf" or "png".
#' @inheritParams verbose logical, whether to print out progress messages and warnings.
#
#' @export
#' @return \code{plotSummaryHeatmap} does not return data but creates image and CSV files. 
#' 
#' 
#' @name plotSummaryHeatmap
#' @aliases plotSummaryHeatmap,EGSEAResults-method
#' @rdname EGSEAResults-methods
#' 
#' @importFrom  RColorBrewer brewer.pal
#' 
#' @examples
#' # Example of plotSummaryHeatmap
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' class(gsa)
#' plotSummaryHeatmap(gsa, gs.label="kegg")
#' 



  setGeneric(name="plotSummaryHeatmap",
          def = function(object, gs.label=1, number = 20, sort.by = NULL,
                  hm.vals = NULL,
                  show.vals = NULL,
                  file.name="sum_heatmap", format = "pdf",
                  verbose=TRUE){
              standardGeneric("plotSummaryHeatmap")
          }
  )
  
  setMethod(f = "plotSummaryHeatmap",
          signature="EGSEAResults",
          definition = function(object, gs.label=1, number = 20, 
                  sort.by = NULL, 
                  hm.vals = NULL,
                  show.vals = NULL,
                  file.name="sum_heatmap", format = "pdf",
                  verbose=TRUE){  
              tryCatch({
              gsc.name = object@gs.annots[[gs.label]]$name
              if (!is.null(sort.by))
                  sort.by = tolower(sort.by)
              else
                  sort.by = object@sort.by
              if (!is.null(hm.vals))
                  hm.vals = tolower(hm.vals)
              else
                  hm.vals = sort.by
              if (length(object@contr.names) > 1){
                  contrast = 0
              }else{
                  contrast = 1
              }   
              if (verbose)
                  cat(paste0("Generating summary heatmap for the collection ",
                          gsc.name, "\n", "sort.by: ",sort.by, ", hm.vals: ",
                          hm.vals, ", show.vals: ", show.vals, "\n"))
      
              t = topSets(object, gs.label, contrast, sort.by, number, TRUE, FALSE)
              hm = matrix(0, length(t), length(object@contr.names))
              if (!is.null(show.vals)){
                  show.vals = tolower(show.vals)
                  cellvals = matrix(0, length(t), length(object@contr.names))                  
                  colnames(cellvals) = paste0(object@contr.names, ".", show.vals)
              }else{
                  cellvals = NULL
              }
              for (i in 1:length(object@contr.names)){
                  hm[,i] = object@results[[gs.label]][["test.results"]][[i]][t, hm.vals]
                  if (!is.null(show.vals))
                      cellvals[,i] = object@results[[gs.label]][["test.results"]][[i]][t, show.vals]
              }
              if (hm.vals %in% c("p.value", "p.adj")){
                  hm = -log10(hm)                  
                  hm[hm == Inf] = max(hm[hm != Inf]) + 50
                  hm[is.na(hm)] = 0 
                  xlab = paste0("-log10(", hm.vals, ")")
              }else
                  xlab = hm.vals
              
              if (length(object@contr.names) > 1){                      
                  colnames(hm) = object@contr.names                  
              }else{
                  hm = cbind(hm, hm)
                  colnames(hm) = c(" ", " ")
                  if (!is.null(show.vals))
                      cellvals = cbind(cellvals, cellvals)                  
              }
              rownames(hm) = t
              title = paste0(gsc.name, " (sorted by ", sort.by, ")")
              suppressWarnings(generateSummaryHeatmaps(
                              hm, hm.vals,                              
                              object@contr.names,
                              title, xlab,
                              file.name, cellvals, format))
              
          }, 
          error = function(e){
              cat(paste0("ERROR: plotSummaryHeatmap(...) encountered an error:\n", e ))
          }) 
          }
  
  )
 
  
#' @title Plot a pathway map for a given KEGG pathway
#' @description \code{plotPathway} generates a visual map for a selected KEGG pathway with
#' the gene fold changes overalid on it.
#' 
#' @inheritParams object EGSEAResults object, the analysis result object from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' @inheritParams gene.set character, the name of the pathway. 
#' See the output of \code{\link{topSets}}.
#' @inheritParams gs.label character or numeric, the index/name of the KEGG pathways collection.
#' See names(object@@gs.annots)[grep("^kegg", names(object@@gs.annots))] for valid values.
#' @inheritParams contrast character or numeric, the index/name of the contrast or 0/"comparison". 
#' See object@@contr.names for valid values.
#' @inheritParams file.name character, the name of the output file without an extension.
#' @inheritParams verbose logical, whether to print out progress messages and warnings.
#' 
#' 
#' @export
#' @return \code{plotPathway} does not return data but creates a file.
#' 
#' 
#' @name plotPathway
#' @aliases plotPathway,EGSEAResults-method
#' @rdname EGSEAResults-methods
#' 
#' @examples
#' # Example of plotPathway
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' class(gsa)
#' plotPathway(gsa, gs.label="kegg", "Asthma")
#' plotPathway(gsa, gs.label="kegg", "Asthma", contrast="comparison", 
#' file.name = "asthma.map.cmp")
#' 
  
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
                                  contrast = object@contr.names[contrast]
                              else
                                  contrast = "comparison"
                          if (verbose)
                              cat(paste0("Generating pathway map for ", gene.set, 
                                              " from the collection \n",
                                              object@gs.annots[[gs.label]]$name, " and for the contrast ",
                                              contrast, "\n"))
                          if (contrast %in% object@contr.names){
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
                                              Use one of the object@contr.names or a numeric value.")
                          }
                      }, 
                      error = function(e){
                          cat(paste0("ERROR: plotPathway(...) encountered an error:\n", e ))
                      })            
          }
  )
  
#' @title Plot a multi-dimensional scaling (MDS) plot for the gene set rankings
#' @description \code{plotMethods} generates a multi-dimensional scaling (MDS) plot
#'  for the gene set rankings of different base GSE methods
#' 
#' @inheritParams object EGSEAResults object, the analysis result object from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' @inheritParams gs.label character or numeric, the index/name of the KEGG pathways collection.
#' See names(object@@gs.annots)[grep("^kegg", names(object@@gs.annots))] for valid values.
#' @inheritParams contrast character or numeric, the index/name of the contrast or 0/"comparison". 
#' See object@@contr.names for valid values.
#' @inheritParams file.name character, the name of the output file without an extension.
#' @inheritParams format character, takes "pdf" or "png".
#' @inheritParams verbose logical, whether to print out progress messages and warnings. 
#' 
#' @export
#' @return \code{plotMethods} does not reutrn data but creates an image file.
#' 
#' 
#' @aliases plotMethods,EGSEAResults-method
#' @rdname EGSEAResults-methods
#'  
#' @examples
#' # Example of plotMethods
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' class(gsa)
#' plotMethods(gsa)
#' 
  
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
              #TODO: color methods based on their null-hypothesis (competitve vs self-contained)
              tryCatch({         
                  if (is.numeric(contrast))
                      if (contrast != 0)
                          contrast = object@contr.names[contrast]
                      else
                          contrast = "comparison"
                  if (tolower(contrast) == "comparison")
                      stop("plotMethods(...) is not supported for the comparison analysis")
                  if (length(object@baseMethods) < 2){
                      stop("plotMethods(...) requires at least two base methods.")
                  }
                  if (verbose)
                      cat(paste0("Generating methods plot for the collection \n",
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
#' @description \code{plotSummary} generates a Summary plot for EGSEA analysis.
#' @details \code{plotSummary} generates a Summmary Plot for an EGSEA analysis.
#' Since the EGSEA "Significance Score" is proportional to the p-value and the 
#' absolute fold changes, it could be useful to highlight gene sets that
#' have high Significance scores. The blue labels on the summary plot indicate 
#' gene sets that do not apear in the top 10 list of gene sets based on the "sort.by" 
#' argument (black labels) yet they appear in the top 5 list of gene sets based on 
#' the EGSEA "Significance Score". If two contrasts are provided, the rank is calculated 
#' based on the "comparison" analysis results and the "Significance Score" is calculated 
#' as the mean. If \code{sort.by = NULL}, the slot \code{sort.by} of the \code{object} 
#' is used to order gene sets.
#' 
#' @inheritParams object EGSEAResults object, the analysis result object from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' @inheritParams gs.label character or numeric, the index/name of the KEGG pathways collection.
#' See names(object@@gs.annots)[grep("^kegg", names(object@@gs.annots))] for valid values.
#' @inheritParams contrast character or numeric, the index/name of the contrast. To compare two
#' contrasts, pass their indexes/names.
#' See object@@contr.names for valid values.
#' @inheritParams file.name character, the name of the heatmap file without an extension.
#' @inheritParams format character, takes "pdf" or "png".
#' @param x.axis character, the x-axis of the summary plot. All the 
#' values accepted by the 
#' \strong{sort.by} parameter can be used. Default x.axis="p.value".
#' @param x.cutoff numeric, cut-off threshold to filter the gene sets of 
#' the summary plots
#' based on the values of the \strong{x.axis}. Default 
#' x.cutoff=NULL.
#' @inheritParams sort.by
#' @param use.names logical, determines whether to display the GeneSet IDs or GeneSet Names.
#' Default is FALSE.
#' @inheritParams verbose logical, whether to print out progress messages and warnings. 
#' 
#' @export
#' @return \code{plotSummary} does not return data but creates an image file. 
#' 
#' 
#' @name plotSummary
#' @aliases plotSummary,EGSEAResults-method
#' @rdname EGSEAResults-methods
#' 
#' @examples
#' # Example of plotSummary
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
                  sort.by = NULL,
                  use.names = FALSE,
                  verbose = TRUE){
              standardGeneric("plotSummary")
          }
  )
  
  setMethod(f = "plotSummary",
          signature="EGSEAResults",
          definition = function(object, gs.label=1, contrast=1, 
                  file.name="summary", format = "pdf",
                  x.axis = "p.adj", x.cutoff = NULL,
                  sort.by = NULL,
                  use.names = FALSE,
                  verbose = TRUE){                
              tryCatch({         
                          if (is.numeric(contrast))                            
                              contrast = object@contr.names[contrast]
                          if (is.null(x.axis))
                              x.axis = object@sum.plot.axis
                          if (length(contrast) == 2){
                              if (verbose)
                                  cat(paste0("Generating Summary plots for the collection \n",
                                                  object@gs.annots[[gs.label]]$name, " and for the comparison ",
                                                  paste(contrast, collapse=" vs "), "\n"))                        
                              ordered.data = object@results[[gs.label]][["comparison"]][["test.results"]]
                              if (!is.null(sort.by))
                                ordered.data = ordered.data[order(ordered.data[, sort.by]), ]
                              plot.data = generatePlotData.comparison(
                                      object@results[[gs.label]][["test.results"]][contrast], 
                                      ordered.data, 
                                      object@gs.annots[[gs.label]], 
                                      x.axis, x.cutoff,
                                      use.names)
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
                              ordered.data = object@results[[gs.label]][["test.results"]][[contrast]]
                              if (!is.null(sort.by))
                                  ordered.data = ordered.data[order(ordered.data[, sort.by]), ]
                              plot.data = generatePlotData(
                                      ordered.data, 
                                      object@gs.annots[[gs.label]], 
                                      x.cutoff, x.axis,
                                      use.names)
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
#' @description \code{plotGOGraph} generates a graph of the top significant GO terms in 
#' a GO term collection, which could be c5 from MSigDB or Gene Ontolog from the GeneSetDB.
#' 
#' @inheritParams object EGSEAResults object, the analysis result object from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' @inheritParams gs.label character or numeric, the index/name of the KEGG pathways collection.
#' See names(object@@gs.annots)[grep("^kegg", names(object@@gs.annots))] for valid values.
#' @inheritParams contrast character or numeric, the index/name of the contrast or 0/"comparison". 
#' See object@@contr.names for valid values.
#' @inheritParams sort.by
#' @param noSig numeric, number of significant GO terms to be displayed. A number larger than 
#' 5 might not work due to the size of the generated graph. 
#' @inheritParams file.name character, the prefix of the output file name without an extension.
#' @inheritParams format
#' @inheritParams verbose logical, whether to print out progress messages and warnings.
 
#' 
#' 
#' @export
#' @return \code{plotGOGraph} does not return data but creates an image file.  
#' 
#' 
#' @name plotGOGraph
#' @aliases plotGOGraph,EGSEAResults-method
#' @rdname EGSEAResults-methods
#' 
#' @examples
#' # Example of plotGOGraph
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' class(gsa)
#' plotGOGraph(gsa, sort.by="avg.rank")
#' 
  
  setGeneric(name="plotGOGraph",
          def = function(object, gs.label="c5", contrast=1, sort.by=NULL,
                  noSig = 5, file.name="c5-top-",
                  format="pdf", verbose=TRUE){
              standardGeneric("plotGOGraph")
          }
  )
  
  setMethod(f = "plotGOGraph",
          signature="EGSEAResults",
          definition = function(object, gs.label="c5", contrast=1, sort.by=NULL,
                  noSig = 5, file.name="GO-top-",
                  format="pdf", verbose=TRUE){            
              tryCatch({               
                          if (is.numeric(gs.label))
                              gs.label = names(object@gs.annots)[gs.label]
                          stopifnot("GOID" %in% colnames(object@gs.annots[[gs.label]]@anno))
                          if (is.numeric(contrast))
                              if (contrast != 0)
                                  contrast = object@contr.names[contrast]
                              else
                                  contrast = "comparison"
                          if (is.null(sort.by))
                              sort.by = object@sort.by
                          if (verbose)
                              cat(paste0("Generating GO Graphs for the collection ",
                                              object@gs.annots[[gs.label]]$name, 
                                              "\n and for the contrast ",
                                              contrast, " based on the ", sort.by, "\n"))
                        
                          if (contrast %in% object@contr.names){
                              results = object@results[[gs.label]][["test.results"]][[contrast]]
                          }else if (tolower(contrast) == "comparison"){
                              results = object@results[[gs.label]][["comparison"]][["test.results"]]
                          }else{
                              stop("Unrecognized contrast value. 
                                              Use one of the object@contr.names or a numeric value.")
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
  
  
#' @title Plot a multi-dimensional scaling (MDS) plot for the gene set rankings
#' @description \code{plotBars} generates a multi-dimensional scaling (MDS) plot
#'  for the gene set rankings of different base GSE methods
#' 
#' @inheritParams object EGSEAResults object, the analysis result object from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' @inheritParams gs.label character or numeric, the index/name of the KEGG pathways collection.
#' See names(object@@gs.annots)[grep("^kegg", names(object@@gs.annots))] for valid values.
#' @inheritParams contrast character or numeric, the index/name of the contrast or 0/"comparison". 
#' See object@@contr.names for valid values.
#' @inheritParams number integer, maximum number of gene sets to list
#' @inheritParams sort.by character
#' @param bar.vals character, determines which EGSEA score values are used to draw the bars.
#' Default is NULL which implies using the sort.by score. 
#' @inheritParams file.name character, the name of the output file without an extension.
#' @inheritParams format character, takes "pdf" or "png".
#' @inheritParams verbose logical, whether to print out progress messages and warnings. 
#' 
#' @export
#' @return \code{plotBars} does not reutrn data but creates an image file.
#' 
#' 
#' @aliases plotBars,EGSEAResults-method
#' @rdname EGSEAResults-methods
#'  
#' @examples
#' # Example of plotBars
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' class(gsa)
#' plotBars(gsa)
#' 
  
  setGeneric(name="plotBars",
          def = function(object, gs.label=1, contrast=1, 
                  number = 20, 
                  sort.by = NULL,
                  bar.vals = "p.adj",
                  file.name="bars_plot", format = "pdf",
                  verbose = TRUE){
              standardGeneric("plotBars")
          }
  )
  
  setMethod(f = "plotBars",
          signature="EGSEAResults",
          definition = function(object, gs.label=1, contrast=1, 
                  number = 20, 
                  sort.by = NULL,
                  bar.vals = "p.adj",
                  file.name="bars_plot", format = "pdf",
                  verbose = TRUE){             
              tryCatch({         
                      if (is.numeric(gs.label))
                          gs.label = names(object@gs.annots)[gs.label]      
                      gsc.name = object@gs.annots[[gs.label]]$name
                      if (is.numeric(contrast))
                          if (contrast != 0)
                              contrast = object@contr.names[contrast]
                          else
                              contrast = "comparison"
    #                          if (tolower(contrast) == "comparison")
    #                              stop("plotBars(...) is not supported for the comparison analysis")
                      if (!is.null(sort.by))
                          sort.by = tolower(sort.by)
                      else
                          sort.by = object@sort.by
                      if (!is.null(bar.vals))
                          bar.vals = tolower(bar.vals)
                      else
                          bar.vals = sort.by  
                      if (verbose)
                          cat(paste0("Generating a bar plot for the collection ",
                                          gsc.name, " \n and the contrast ",
                                          contrast, "\n"))                          
                      t = topSets(object, gs.label, contrast, sort.by, number, FALSE, FALSE)
                      t = t[nrow(t):1, ]
                      vals = t[, bar.vals]
                      if (bar.vals %in% c("p.value", "p.adj")){
                          vals = -log10(vals)
                          vals[vals == Inf] = max(vals[vals != Inf]) + 50
                          vals[is.na(vals)] = 0 
                          xlab = paste0("-log10(", bar.vals, ")")
                      }else
                          xlab = bar.vals
                      col = rep("#67A9CF", length(vals))
                      col[t[, "direction"] == "Up"] = "#EF8A62"
                      col[t[, "direction"] == "Neutral"] = "#87669A"
                      sel.genes.sets = length(vals)
                      if(sel.genes.sets <= 20){
                          cr = 0.7
                      }else if(sel.genes.sets < 40){
                          cr = 0.65
                      }else if(sel.genes.sets < 70){
                          cr = 0.35
                      }else if(sel.genes.sets < 100){
                          cr = 0.35
                      }else 
                          cr = 0.15
                      title = paste0(gsc.name, " (sorted by ", sort.by, ")")
                      t1 = rownames(t)                                                    
                      for (i in 1:length(t1)){                      
                          if (nchar(t1[i]) > 25)
                              t1[i] = paste0(substr(t1[i], 1, 25), " ...")
                      }                         
                      if (is.null(format) || tolower(format) == "pdf"){
                          pdf(paste0(file.name, ".pdf"))
                          par(las=1)
                          par(mar=c(5,12,4,2))
                          par(cex.main = 0.8)
                          if (abs(max(vals)) != Inf)
                            xlim = c(0, max(vals) * 1.5)
                          else
                            xlim = NULL  
                          barplot2(vals, horiz = TRUE, 
                                  names.arg = t1,
                                  cex.names = cr,
                                  col = col,
                                  xlab = xlab, 
                                  xlim = xlim, 
                                  main = title)   
                          if (bar.vals %in% c("p.value", "p.adj")){
                              abline(v = -log10(0.05), col = "red", lwd = 1, lty=2)
                          }
                          legend("topright", xpd=TRUE,   
                                  legend = c("Up-regulated", 
                                          "Neutral", 
                                          "Down-regulated"),
                                  border = "#FFFFFF",
                                  fill = "#FFFFFF",
                                  col = c("#EF8A62", "#87669A", "#67A9CF"), 
                                  title = "Regulation Direction",                             
                                  lty= 1,             
                                  lwd = 5,           
                                  cex=.7
                          )
                          dev.off()
                      }
                      if (is.null(format) || tolower(format) == "png"){
                          png(paste0(file.name, ".png"))
                          par(las=1)
                          par(mar=c(5,12,4,2))
                          par(cex.main = 0.8)
                          barplot2(vals, horiz = TRUE, 
                                  names.arg = t1,
                                  cex.names = cr,
                                  col = col,
                                  xlab = xlab, 
                                  main = title)   
                          if (bar.vals %in% c("p.value", "p.adj")){
                              abline(v = -log10(0.05), col = "red", lwd = 1, lty=2)
                          }
    #                              legend("topright", xpd=TRUE,   
    #                                      legend = c("Up-regulated", "Down-regulated"),
    #                                      border = "#FFFFFF",
    #                                      fill = "#FFFFFF",
    #                                      col = c("#EF8A62", "#67A9CF"), 
    #                                      title = "Regulation Direction",                             
    #                                      lty= 1,             
    #                                      lwd = 5,           
    #                                      cex=.7
    #                              )
                          dev.off()
                      }
                  }, 
                  error = function(e){
                      cat(paste0("ERROR: plotBars(...) encountered an error:\n", e ))
                  }
              )
              
          }
  )
  
  
#' @title Display the information of a given gene set name
#' @description \code{showSetByname} shows the details of a given gene set indicated by name.
#' 
#' @inheritParams object EGSEAResults 
#' @inheritParams gs.label
#' @param set.name character, a vector of gene set names as they appear in \code{\link{topSets}}.
#' 
#' @export
#' @return \code{showSetByName} does not return data
#' 
#' 
#' @name showSetByName
#' @aliases showSetByName,EGSEAResults-method
#' @rdname EGSEAResults-methods
#' 
#' @examples
#' # Example of showSetByName
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' class(gsa)
#' showSetByName(gsa, "kegg", "Asthma")
#' 
  
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
#' @description \code{showSetByID} shows the details of a given gene set indicated by ID.
#' 
#' @inheritParams object EGSEAResults
#' @inheritParams gs.label
#' @param id character, a vector of gene set IDs as they appears in the 
#' \code{\link{plotSummary}}.
#' 
#' @export
#' @return  \code{showSetByID} does not return data.
#' 
#' 
#' @name showSetByID
#' @aliases showSetByID,EGSEAResults-method
#' @rdname EGSEAResults-methods
#'  
#' @examples
#' # Example of showSetByID
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' class(gsa)
#' showSetByID(gsa, "kegg", "hsa04060")
#' 
  
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
  

#' @title The gene set enrichment scores per sample
#' @description \code{getSetScores} returns a dataframe of the gene set enrichment scores 
#' per sample. This can be only calculated using specific base methods, namely, "ssgsea". 
#' 
#' @inheritParams object EGSEAResults
#' @inheritParams gs.label
#' 
#' @export
#' @return  \code{getSetScores} returnsa a dataframe where rows are gene sets and 
#' columns are samples.
#' 
#' 
#' @name getSetScores
#' @aliases getSetScores,EGSEAResults-method
#' @rdname EGSEAResults-methods
#'  
#' @examples
#' # Example of getSetScores
#' library(EGSEAdata)
#' data(il13.gsa)
#' gsa = il13.gsa
#' class(gsa)
#' head(getSetScores(gsa, "kegg"))
#' 
  
  setGeneric(name="getSetScores",
          def = function(object, gs.label = 1){ # , method = "ssgsea"
              standardGeneric("getSetScores")
          }
  )
# @param method character, a method that calculates gene set enrichment score per sample.
# Supported method is "ssgsea".
  setMethod(f = "getSetScores",
          signature(object = "EGSEAResults"),
          definition = function(object, gs.label = 1){ # , method = "ssgsea"
              method = "ssgsea"
              if (is.numeric(gs.label))
                  stopifnot(gs.label > 0 && gs.label <= length(object@gs.annots))
              else
                  stopifnot(gs.label %in% names(object@gs.annots))
              if ("set.scores" %in% names(object@results[[gs.label]])){
                  method = tolower(method)
                  stopifnot(method %in% object@baseMethods)
                  return(object@results[[gs.label]][["set.scores"]][[method]])
              }else{
                  cat("set.scores are not available for this EGSEAResults object.\n")
                  cat("Try to re-run EGSEA analysis with cal.set.scores = TRUE \n")
                  cat("and use one of the set scoring methods, i.e., ssgsea.")
                  return(NULL)
              }
          }
  
  )
  
  