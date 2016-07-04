#Ensemble of Gene Set Enrichment Analyses 
# 
# Author: Monther Alhamdoosh, E:m.hamdoosh@gmail.com

###############################################################################


#' The GSCollectionIndex class
#'
#' This class contains an indexed gene set collection
#'
#' @slot original list, the original gene sets
#' @slot idx list, the gene set indexes
#' @slot anno data.frame, the annotations of the gene sets
#' @slot featureIDs character, vector of the original Entrez IDs that are 
#' used in the indexing procedure 
#' @slot species character, the species name 
#' @slot name character, the name of the gene set collection
#' @slot label character, a label to distnguish this collection 
#' from other collections 
#' 
#' 
#' @name GSCollectionIndex 
#' @rdname GSCollectionIndex-methods
#' @aliases GSCollectionIndex-class
#' @exportClass GSCollectionIndex

GSCollectionIndex <- setClass(
        "GSCollectionIndex",            
        slots = c(original = "list",
                idx = "list",
                anno = "data.frame",                
                featureIDs = "character",
                species = "character",
                name = "character",
                label = "character"),                      
        prototype = list(original = list(),
                idx = list(),
                anno = data.frame(),                
                featureIDs = c(),
                species = "",
                name = "",
                label = "")           
)

#' @title Extract a slot from an object of class GSCollectionIndex
#' @description Extract a slot from an object of class GSCollectionIndex
#' 
#' @param x GSCollectionIndex, the indexed gene set collection generated
#' from \code{\link{buildIdx}}, \code{\link{buildMSigDBIdx}}, 
#' \code{\link{buildKEGGIdx}},
#' \code{\link{buildGeneSetDBIdx}}, and \code{\link{buildCustomIdx}}. 
#' @param name character, the slot name
#' 
#' @export
#' @return the selected slot.
#' 
#' 
#' @aliases $,GSCollectionIndex-method
#' @rdname GSCollectionIndex-methods
#' 
#' @examples
#' library(EGSEAdata)
#' data(il13.data)
#' v = il13.data$voom
#' gs.annots = buildIdx(entrezIDs=rownames(v$E), species="human", 
#' 			msigdb.gsets="none", 
#'          kegg.updated=FALSE, kegg.exclude = c("Metabolism"))
#' 
#' print(gs.annots[[1]]$name)

setMethod("$", "GSCollectionIndex",
        function(x, name){           
            slot(x, name)
        }
)

setGeneric(name="selectGeneSets",
        def = function(object, gs.names=NULL, min.size=1){
            standardGeneric("selectGeneSets")
        }
)

setMethod(f = "selectGeneSets",
        signature(object = "GSCollectionIndex"),
        definition = function(object, gs.names=NULL, min.size=1){
            if (length(object@idx) == 0)
                return(object)
            gs.annot.top = GSCollectionIndex()
            if (is.null(gs.names)){
                gs.names = names(object@idx[sapply(object@idx, function(x) 
                                            length(x)) >= min.size])    
            }   
            else{
                gs.names = gs.names[sapply(object@idx[gs.names], function(x) 
                                    length(x)) >= min.size]
            }
            sel = match(gs.names, object@anno[, "GeneSet"])
            gs.annot.top@original = object@original[sel]
            gs.annot.top@idx = object@idx[sel]
            gs.annot.top@anno = object@anno[sel,]
            gs.annot.top@label = object@label
            gs.annot.top@featureIDs = object@featureIDs
            gs.annot.top@species = object@species
            gs.annot.top@name = object@name
            return(gs.annot.top)          
        }

)


#' @title Print a summary of the gene set collection
#' @description A brief summary of the gene set collection
#' 
#' @inheritParams object GSCollectionIndex, the indexed gene set collection generated
#' from \code{\link{buildIdx}}, \code{\link{buildMSigDBIdx}}, 
#' \code{\link{buildKEGGIdx}},
#' \code{\link{buildGeneSetDBIdx}}, and \code{\link{buildCustomIdx}}. 
#' 
#' @export
#' @return summary of the collection.
#' 
#' 
#' @aliases summary,GSCollectionIndex-method
#' @rdname GSCollectionIndex-methods
#'  
#' @examples
#' library(EGSEAdata)
#' data(il13.data)
#' v = il13.data$voom
#' gs.annots = buildIdx(entrezIDs=rownames(v$E), species="human", 
#' 			msigdb.gsets="none", 
#'          kegg.updated=FALSE, kegg.exclude = c("Metabolism"))
#' 
#' summary(gs.annots[[1]])

setMethod(f = "summary",
        signature(object = "GSCollectionIndex"),
        definition = function(object){
            cat(paste0(object@name, " (", object@label , "): ",
                                  length(object@idx), " gene sets"))
        }

)



#' @title Display the content of the gene set collection
#' @description Print the details of the gene set collection
#' 
#' @inheritParams object GSCollectionIndex, the indexed gene set collection generated
#' from \code{\link{buildIdx}}, \code{\link{buildMSigDBIdx}}, 
#' \code{\link{buildKEGGIdx}},
#' \code{\link{buildGeneSetDBIdx}}, and \code{\link{buildCustomIdx}}. 
#' 
#' @export
#' @return details of the collection.
#' 
#' 
#' @aliases show,GSCollectionIndex-method
#' @rdname GSCollectionIndex-methods
#'  
#' @examples
#' library(EGSEAdata)
#' data(il13.data)
#' v = il13.data$voom
#' gs.annots = buildIdx(entrezIDs=rownames(v$E), species="human", 
#' 			msigdb.gsets="none", 
#'          kegg.updated=FALSE, kegg.exclude = c("Metabolism"))
#' 
#' show(gs.annots[[1]])

setMethod(f = "show",
        signature(object = "GSCollectionIndex"),
        definition = function(object){
            cat("An object of class \"GSCollectionIndex\"\n")
            cat(paste0("\tNumber of gene sets: ", length(object@original), "\n"))
            if (nrow(object@anno) != 0)
                cat(paste0("\tAnnotation columns: ", 
                                paste(colnames(object@anno), collapse=", "), "\n"))
            else
                cat(paste0("\tAnnotation not provided", "\n"))
            cat(paste0("\tTotal number of indexing genes: ", 
                            length(object@featureIDs), "\n"))
            cat(paste0("\tSpecies: ", object@species, "\n"))
            cat(paste0("\tCollection name: ", object@name, "\n"))
            cat(paste0("\tCollection uniqe label: ", object@label, "\n"))            
        }

)



#' @title Display the information of a given gene set name
#' @description Print the details of a given gene set name
#' 
#' @param object GSCollectionIndex, the indexed gene set collection generated
#' from \code{\link{buildIdx}}, \code{\link{buildMSigDBIdx}}, 
#' \code{\link{buildKEGGIdx}},
#' \code{\link{buildGeneSetDBIdx}}, and \code{\link{buildCustomIdx}}. 
#' @param set.name character, a vector of gene set names as they appear in \code{\link{topSets}}.
#' 
#' @export
#' @return a list of the details of a gene set
#' 
#' 
#' @name getSetByName
#' @aliases getSetByName,GSCollectionIndex-method
#' @rdname GSCollectionIndex-methods
#' 
#' @examples
#' library(EGSEAdata)
#' data(il13.data)
#' v = il13.data$voom
#' gs.annots = buildIdx(entrezIDs=rownames(v$E), species="human", 
#' 			msigdb.gsets="none", 
#'          kegg.updated=FALSE, kegg.exclude = c("Metabolism"))
#' 
#' getSetByName(gs.annots[[1]], "Asthma")

setGeneric(name="getSetByName",
        def = function(object, set.name){
            standardGeneric("getSetByName")
        }
)

setMethod(f = "getSetByName",
        signature(object = "GSCollectionIndex"),
        definition = function(object, set.name){
            stopifnot(set.name %in% names(object@original))
            cols = colnames(object@anno)
            gset = list()
            for (x in set.name){
                gset[[x]] = list()
                for (col in cols){
                    cat(paste0(col, ": ", object@anno[x, col], "\n"))
                    gset[[col]] = object@anno[x, col]
                }
                cat("\n")                
            }
            return(gset)
        }

)



#' @title Display the information of a given gene set ID
#' @description Print the details of a given gene set ID
#' 
#' @inheritParams object GSCollectionIndex, the indexed gene set collection generated
#' from \code{\link{buildIdx}}, \code{\link{buildMSigDBIdx}}, 
#' \code{\link{buildKEGGIdx}},
#' \code{\link{buildGeneSetDBIdx}}, and \code{\link{buildCustomIdx}}. 
#' @param id character, a vector of gene set IDs as they appears in the 
#' \code{\link{plotSummary}}.
#' 
#' @export
#' @return  a list of the details of a gene set
#' 
#' 
#' @name getSetByID
#' @aliases getSetByID,GSCollectionIndex-method
#' @rdname GSCollectionIndex-methods
#'  
#' @examples
#' library(EGSEAdata)
#' data(il13.data)
#' v = il13.data$voom
#' gs.annots = buildIdx(entrezIDs=rownames(v$E), species="human", 
#' 			msigdb.gsets="none", 
#'          kegg.updated=FALSE, kegg.exclude = c("Metabolism"))
#' 
#' getSetByID(gs.annots[[1]], "hsa04060")

setGeneric(name="getSetByID",
        def = function(object, id){
            standardGeneric("getSetByID")
        }
)

setMethod(f = "getSetByID",
        signature(object = "GSCollectionIndex"),
        definition = function(object, id){
            stopifnot("ID" %in% colnames(object@anno))
            stopifnot(id %in% object@anno$ID)            
            name = which(object@anno$ID %in% id)
            cols = colnames(object@anno)
            gset = list()
            for (x in name){
                gset[[x]] = list()
                for (col in cols){
                    cat(paste0(col, ": ", object@anno[x, col], "\n"))
                    gset[[col]] = object@anno[x, col]
                }
                cat("\n")
            }
            return(gset)
        }

)