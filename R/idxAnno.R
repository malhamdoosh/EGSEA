#Ensemble of Gene Set Enrichment Analyses 
# 
# Author: Monther Alhamdoosh, E:m.hamdoosh@gmail.com
###############################################################################


species.fullToShort = list()
species.fullToShort[["homo sapiens"]] = "human"
species.fullToShort[["mus musculus"]] = "mouse"
species.fullToShort[["rattus norvegicus"]] = "rat"
# "Homo sapiens", "Mus musculus", "Rattus norvegicus","Danio rerio","Macaca 
# mulatta"

msigdb.gs.names = list(c1="c1 Positional Gene Sets", c2="c2 Curated Gene Sets", 
c3="c3 Motif Gene Sets", 
        c4="c4 Computational Gene Sets", c5="c5 GO Gene Sets", 
        c6="c6 Oncogenic Signatures", c7="c7 Immunologic Signatures",
        h="h Hallmark Signatures")
genesetdb.gs.labels = list("GeneSetDB Drug/Chemical"="gsdbdrug",   
 "GeneSetDB Disease/Phenotype" = "gsdbdis",
 "GeneSetDB Gene Ontology" = "gsdbgo",
 "GeneSetDB Pathway" = "gsdbpath",
 "GeneSetDB Gene Regulation" = "gsdbreg")
human.names = c("human", "homo sapiens", "hs")
mouse.names = c("mouse", "mus musculus", "mm")
rat.names = c("rat", "rattus norvegicus" , "rn")
loadKeggData <- function(){
    data("kegg.pathways", package="EGSEAdata")
    if (is.null(kegg.pathways))
        stop("Failed to load the KEGG pathway data.")
    return (kegg.pathways)
}

#' @title Gene Set Collection Index from the KEGG Database
#' 
#' @description It prepares the KEGG pathway collection to be used for the 
#' EGSEA analysis. 
#' 
#' @details It indexes the KEGG pathway gene sets and loads gene set annotation.
#'   
#' @param entrezIDs character, a vector that stores the Entrez Gene IDs tagged 
#' in your dataset. 
#' The order of the Entrez Gene IDs should match those of the count/expression 
#' matrix row names.
#' @param species character, determine the organism of selected gene sets: 
#' "human", "mouse" or "rat".
#' @param min.size integer, the minium number of genes required in a testing 
#' gene set
#' @param updated logical, set to TRUE if you want to download the most recent 
#' KEGG pathways.
#' @param rdata.dir character, directory from which the KEGG pathway 
#' collections are loaded. 
#' If NULL, EGSEA tries to load the data from \pkg{EGSEAdata}. 
#'
#' @return indexed gene set annotation that can be used with other functions in 
#' the package.
#' Each annotation is a list of seven elements: \code{original} stores the 
#' original gene sets, 
#' \code{idx} stores the indexed gene sets,  \code{anno} that stores detailed 
#' annotation for each 
#' gene set, \code{label} a unique id that identifies the collection of gene 
#' sets, 
#' \code{featureIDs} stores the entrezIDs used in building the annotation, 
#' \code{species}
#'  stores that organism name of gene sets and \code{name}  stores the 
#' collection name 
#' to be used in the EGSEA report.
#' 
#' @import EGSEAdata
#' @export 
#' @examples 
#' library(EGSEAdata)
#' data(il13.data)
#' v = il13.data$voom
#' gs.annots = buildKEGGIdxEZID(entrezIDs=rownames(v$E), species="human")
#' 

buildKEGGIdxEZID <- function(entrezIDs, species = "human", min.size=1, 
updated=FALSE, 
        rdata.dir=NULL) {   
    species = tolower(species)
    if (species %in% human.names){
        species = "Homo sapiens"
    }else if (species %in% mouse.names){
        species = "Mus musculus"
    }
    else if (species %in% rat.names){
        species = "Rattus norvegicus"
    }
    else{
        stop("Unrecognized species.")
    }
    entrezIDs = as.character(entrezIDs)
    kegg.file = paste0(rdata.dir, "/kegg.pathways.rda")
    kegg = NULL
    if (!is.null(rdata.dir) && file.exists(kegg.file)){
        print(paste0("KEGG pathways have been loaded from ", kegg.file))
        load(kegg.file)     
        kegg = kegg.pathways[[species.fullToShort[[tolower(species)]]]]
    }else {
        print("Building KEGG pathways annotation object ... ")      
        if (!updated){
            kegg.pathways = loadKeggData()
            kegg = 
kegg.pathways[[species.fullToShort[[tolower(species)]]]]
        }else{
            kegg = tryCatch({                   
    
                        kegg.gsets(species = 
species.fullToShort[[tolower(species)]], id.type = "kegg")          
        
                    },
                    error = function(e){    
                        warning("KEGG pathways have not 
been updated successfully.")
                        #kegg.pathways = NULL
                        kegg.pathways = loadKeggData()  
                                        
                        
return(kegg.pathways[[species.fullToShort[[tolower(species)]]]])
                    })
        }
    }
    if (is.null(kegg))
        stop("Failed to load the KEGG pathway data.")
    gsets = kegg$kg.sets        
    gsets.ez = gsets
    gsets = ids2indices(gsets.ez, entrezIDs, remove.empty=TRUE) 
    gsets.ez = gsets.ez[names(gsets)]   
    gsets.ids = sapply(names(gsets.ez), function (x) as.character(substr(x, 
1,8))) 
    gsets.names = sapply(names(gsets.ez), function (x) 
as.character(substr(x, 10, nchar(x))))
    
    tmp = rep(NA, length(kegg$kg.sets))
    tmp[kegg$sig.idx] = "Signaling"
    tmp[kegg$met.idx] = "Metabolism"
    tmp[kegg$dise.idx] = "Disease"
    anno = data.frame(ID=gsets.ids, GeneSet=gsets.names,
            NumGenes=paste0(sapply(gsets, length), "/", 
sapply(gsets.ez, length)), 
            Type=tmp[match(names(gsets), names(kegg$kg.sets))])
    
    rownames(anno) = gsets.names
    names(gsets) = gsets.names
    names(gsets.ez) = gsets.names
    
    gs.annot = list(original=gsets.ez, idx=gsets, anno=anno, label="kegg", 
featureIDs=entrezIDs,
            species=species, name="KEGG Pathways")
    gs.annot = getGsetAnnot(gs.annot=gs.annot, min.size=min.size)   
    if (length(gs.annot$idx) == 0)
        print("KEGG pathway collection is empty.")
    
    return(gs.annot)
}


#' @title Gene Set Collection Indexes from the MSigDB Database
#' 
#' @description It prepares the MSigDB gene set collections to be used for the 
#' EGSEA analysis. 
#' 
#' @details It indexes the MSigDB gene sets and loads gene set annotation.
#'   
#' @param entrezIDs character, a vector that stores the Entrez Gene IDs tagged 
#' in your dataset. 
#' The order of the Entrez Gene IDs should match those of the count/expression 
#' matrix row names.
#' @param geneSets character, a vector determines which gene set collections 
#' should be used.
#' It can take values from this list: "h", "c1", "c2", "c3", "c4", "c5", 
#' "c6","c7". "h" and "c1"
#' are human specific. If NULL, all available gene set collections are loaded.  
#' @param species character, determine the organism of selected gene sets: 
#' "human", "mouse" or "rat". 
#' @param min.size integer, the minium number of genes required in a testing 
#' gene set
#' @param rdata.dir character, directory from which the MSigDB collections are 
#' loaded. 
#' If NULL, EGSEA tries to load the data from \pkg{EGSEAdata}. 
#'
#'@return indexed gene set annotation that can be used with other functions in 
#' the package.
#' Each annotation is a list of seven elements: \code{original} stores the 
#' original gene sets, 
#' \code{idx} stores the indexed gene sets,  \code{anno} that stores detailed 
#' annotation for each 
#' gene set, \code{label} a unique id that identifies the collection of gene 
#' sets, 
#' \code{featureIDs} stores the entrezIDs used in building the annotation, 
#' \code{species}
#'  stores that organism name of gene sets and \code{name}  stores the 
#' collection name 
#' to be used in the EGSEA report.
#' 
#' @import EGSEAdata
#' @export 
#' @examples 
#' library(EGSEAdata)
#' data(il13.data)
#' v = il13.data$voom
#' gs.annots = buildMSigDBIdxEZID(entrezIDs=rownames(v$E), geneSets=c("h", 
#' "c2"), species="human")
#' names(gs.annots)


buildMSigDBIdxEZID <- function(entrezIDs, geneSets=NULL, species="Homo 
sapiens", 
        min.size=1, rdata.dir=NULL){
    species = tolower(species)
    if (species %in% human.names){
        species = "Homo sapiens"
        if (is.null(geneSets))
            geneSets = c("h", "c1", "c2", "c3", "c4", "c5", 
"c6","c7")
    }else if (species %in% mouse.names){
        species = "Mus musculus"
        if (is.null(geneSets))
            geneSets = c("h", "c2", "c3", "c4", "c5", "c6","c7")
    }
    else{
        stop("Unrecognized species. MSigDB supports human and mouse 
only.")
    }
    geneSets = tolower(geneSets)
    entrezIDs = as.character(entrezIDs)
    gs.annots = list() # Gene set annotation indexes
    gs.annots.names=c() 
    # Read all gene set collections
    print("Reading Broad Gene Sets ... ")
    # SymbolIdentifier(), EntrezIdentifier()
    gsc.all = NULL
    msigdb.file = paste0(rdata.dir, "/msigdb_v5.rda")
    if (!is.null(rdata.dir) && file.exists(msigdb.file))
        load(msigdb.file)
    else
        data("msigdb_v5", package="EGSEAdata")
    gsc.all = gsc.all5              
    if (is.null(gsc.all))
        stop("Failed to load the MSigDB gene set collection data")
    gsc.all = gsc.all[names(gsc.all) == "GENESET"]          
    names(gsc.all) = sapply(gsc.all, function(x) x["STANDARD_NAME"])
    organisms = sapply(gsc.all, function(x) x["ORGANISM"])
    if (species == "Homo sapiens")
        gsc.all = gsc.all[organisms == "Homo sapiens"] #  "Homo 
#sapiens","Mus musculus", "Rattus norvegicus","Danio rerio","Macaca mulatta"
    types = tolower(sapply(gsc.all, function(x) x["CATEGORY_CODE"])) 
#character vector c1 c2 c3 c4 or c5
    ## process each gene set collection
    for (geneSet in geneSets){      
        gs.annots.names = c(gs.annots.names, geneSet)           
                        
        gs.annot = list()
        gsc.small = gsc.all[types == geneSet] 
        if (species == "Mus musculus"){
            if (geneSet %in% c("h", "c2", "c3", "c4", "c5", 
"c6","c7")){
                ver = ifelse (geneSet == "c5", "4", "5")    
            
                msigdb.file = paste0(rdata.dir, 
"/msigdb.mouse_", geneSet,"_v", ver, "5.rdata")
                if (!is.null(rdata.dir) && 
file.exists(msigdb.file))
                    load(msigdb.file)
                else
                    data(list=paste0("msigdb.mouse_", 
geneSet,"_v", ver), package="EGSEAdata")
                geneSet = ifelse(geneSet == "h", "H", geneSet)
                gs.annot$original =get(paste0("Mm.", geneSet))
            }
            else{
                warning(paste0("Unsupported gene set for Mus 
musculus ... ", geneSet))
                next
            }
        }
        else if (species == "Homo sapiens")
            gs.annot$original = sapply(gsc.small, function (x) 
strsplit(x["MEMBERS_EZID"], ",")[[1]])      
        print(paste0("Created the gs.annot$original for ", geneSet, " 
..."))
        gs.annot$idx = ids2indices(gs.annot$original, entrezIDs, 
remove.empty=TRUE) # pathways
        if (length(gs.annot$idx) == 0){
            print(paste0("None of the genes in ", geneSet, "are 
mapped to your gene IDs"))          
        }
        else{
            print(paste0("Created the gs.annot$idx for ", geneSet, 
" ..."))            
            gsc.small = gsc.small[names(gs.annot$idx)] # remove 
#empty sets
            gs.annot$original = 
gs.annot$original[names(gs.annot$idx)]
            
            tmp = rep(NA, length(gs.annot$idx)) 
            gs.annot$anno = data.frame(ID=tmp, 
GeneSet=names(gs.annot$idx), BroadUrl=tmp, Description=tmp,
                    PubMedID=tmp, NumGenes=rep(0, 
length(gs.annot$idx)), Contributor=tmp)
            gs.annot$anno[, "ID"] = sapply(gsc.small, function(x) 
x["SYSTEMATIC_NAME"])
            gs.annot$anno[, "BroadUrl"] = 
paste0("http://www.broadinstitute.org/gsea/msigdb/cards/", names(gs.annot$idx), 
                    ".html")    
            gs.annot$anno[, "Description"] = sapply(gsc.small, 
function(x) x["DESCRIPTION_BRIEF"])
            gs.annot$anno[, "PubMedID"] = sapply(gsc.small, 
function(x) x["PMID"])
            gs.annot$anno[, "NumGenes"] = 
paste0(sapply(gs.annot$idx, length), '/',sapply(gs.annot$original, length)) 
    
            gs.annot$anno[, "Contributor"] = sapply(gsc.small, 
function(x) x["CONTRIBUTOR"])   # KEGG      
            rownames(gs.annot$anno) = names(gs.annot$original)
            print(paste0("Created the gs.annot$anno for ", geneSet, 
" ..."))            
        }
        
        gs.annot$label = geneSet    
        gs.annot$featureIDs = entrezIDs
        gs.annot$species = species
        gs.annot$name = msigdb.gs.names[[gs.annot$label]]       
        if (geneSet == "c5"){
#           print(colnames(gs.annot$anno))
            go.terms = gsub(".*(GO:[0-9]+).*$","\\1", 
gs.annot$anno[, "Description"]) 
            xx = GOTERM
            ontology = character(0)
            sel.gsets = character(0)
            goID = character(0)
            for (i in 1:length(go.terms)){              
                temp = xx[[go.terms[i]]]
                if (!is.null(temp)){
                    ontology = c(ontology, Ontology(temp))
                    sel.gsets = c (sel.gsets, 
as.character(gs.annot$anno[i, "GeneSet"]))
                    goID = c(goID, go.terms[i])
                }
            }
            gs.annot = getGsetAnnot(gs.annot, sel.gsets)
            gs.annot$anno = cbind(gs.annot$anno, Ontology=ontology, 
GOID=goID)
            gs.annot$anno = droplevels(gs.annot$anno)       
    
        }
        gs.annot = getGsetAnnot(gs.annot=gs.annot, min.size=min.size)   
        if (length(gs.annot$idx) == 0)
            print(paste0("MSigDB ", gs.annot$label, " gene set 
collection is empty."))
        gs.annots[[length(gs.annots.names)]] = gs.annot
    }   
    names(gs.annots) = gs.annots.names
    return(gs.annots)
    
}


#' @title Gene Set Collection Indexes from the GeneSetDB Database
#' 
#' @description It prepares the GeneSetDB gene set collections to be used for 
#' the EGSEA analysis. 
#' 
#' @details It indexes the GeneSetDB gene sets and loads gene set annotation.
#'   
#' @param entrezIDs character, a vector that stores the Entrez Gene IDs tagged 
#' in your dataset. 
#' The order of the Entrez Gene IDs should match those of the count/expression 
#' matrix row names.
#' @param species character, determine the organism of selected gene sets: 
#' "human", "mouse" or "rat".
#' @param by.category logical, whether to group the gene sets into collections 
#' or not. Default TRUE.
#' @param min.size integer, the minium number of genes required in a testing 
#' gene set
#' @param rdata.dir character, directory from which the GeneSetDB collections 
#' are loaded. 
#' If NULL, EGSEA tries to load the data from \pkg{EGSEAdata}. 
#'
#' @return indexed gene set annotation that can be used with other functions in 
#' the package.
#' Each annotation is a list of seven elements: \code{original} stores the 
#' original gene sets, 
#' \code{idx} stores the indexed gene sets,  \code{anno} that stores detailed 
#' annotation for each 
#' gene set, \code{label} a unique id that identifies the collection of gene 
#' sets, 
#' \code{featureIDs} stores the entrezIDs used in building the annotation, 
#' \code{species}
#'  stores that organism name of gene sets and \code{name}  stores the 
#' collection name 
#' to be used in the EGSEA report.
#' 
#' @import EGSEAdata
#' @export 
#' @examples 
#' library(EGSEAdata)
#' data(il13.data)
#' v = il13.data$voom
#' gs.annots = buildGeneSetDBIdxEZID(entrezIDs=rownames(v$E), species="human")
#' names(gs.annots)
#' 
#' 


buildGeneSetDBIdxEZID <- function(entrezIDs, species, by.category=TRUE, 
min.size=1, rdata.dir=NULL){
    
    species = tolower(species)  
    if (species %in% human.names){
        species = "Homo sapiens"
    }else if (species %in% mouse.names){
        species = "Mus musculus"
    }
    else if (species %in% rat.names){
        species = "Rattus norvegicus"
    }
    else{
        stop("Unrecognized species.")
    }
    gsdb.file = paste0(rdata.dir, '/genesetdb.', 
species.fullToShort[[tolower(species)]], ".rdata")
    if (!is.null(rdata.dir) && file.exists(gsdb.file))
        load(gsdb.file)
    else
        data(list=paste0('genesetdb.', 
species.fullToShort[[tolower(species)]]), 
                package="EGSEAdata")
    if (is.null(gsetdb.all)){
        stop("Failed to load the GeneSetDB collection data")
    }
    gs.annot = buildCustomIdxEZID(entrezIDs = entrezIDs, gsets = 
gsetdb.all$original,
            anno = gsetdb.all$anno, label = "gsDB", name="GeneSetDB 
Gene Sets", 
            species = species)
    if (by.category){
        categories = unique(as.character(gs.annot$anno$Category))
        gs.annots = list()
        for (cat in categories){
            label = genesetdb.gs.labels[[cat]]
            tmp.gs.annot = getGsetAnnot(gs.annot = gs.annot, 
                    gs.names = 
as.character(gs.annot$anno[gs.annot$anno[,"Category"] == cat,"GeneSet"]),  
                    min.size=min.size)
            if (length(tmp.gs.annot$idx) == 0){
                print(paste0("GeneSetDB ", label," gene set 
collection is empty."))             
            }
            gs.annots[[label]] = tmp.gs.annot
            gs.annots[[label]]$label = label
            gs.annots[[label]]$name = cat
            gs.annots[[label]]$anno = gs.annots[[label]]$anno[, 
-c(match("Category", names(gs.annot$anno)))]
        }
    }else{
        gs.annot = getGsetAnnot(gs.annot=gs.annot, min.size=min.size)   
        if (length(gs.annot$idx) == 0)
            print("GeneSetDB gene set collection is empty.")
        gs.annots = gs.annot
    }
    return(gs.annots)
}

getGsetAnnot <- function(gs.annot, gs.names=NULL, min.size=1){
#   gs.annot = list(original=gsets.sym, idx=gsets, anno=anno, label="kegg", 
#egIds=egIds, featureIDs=symbols, species=species)
    if (length(gs.annot$idx) == 0)
        return(gs.annot)
    gs.annot.top = list()
    if (is.null(gs.names)){
        gs.names = names(gs.annot$idx[sapply(gs.annot$idx, function(x) 
length(x)) >= min.size])    
    }   
    else{
        gs.names = gs.names[sapply(gs.annot$idx[gs.names], function(x) 
length(x)) >= min.size]
    }
    sel = match(gs.names, gs.annot$anno[, "GeneSet"])
    gs.annot.top$original = gs.annot$original[sel]
    gs.annot.top$idx = gs.annot$idx[sel]
    gs.annot.top$anno = gs.annot$anno[sel,]
    gs.annot.top$label = gs.annot$label
    gs.annot.top$featureIDs = gs.annot$featureIDs
    gs.annot.top$species = gs.annot$species
    gs.annot.top$name = gs.annot$name
    return(gs.annot.top)    
}


#' @title Custom Gene Set Collection Index 
#' 
#' @description It creates gene set collections from a given list of gene sets 
#' to be used 
#' for the EGSEA analysis. 
#' 
#' @details It indexes newly created gene sets and attach gene set annotation 
#' if provided.
#'   
#' @param entrezIDs character, a vector that stores the Entrez Gene IDs tagged 
#' in your dataset. 
#' The order of the Entrez Gene IDs should match those of the count/expression 
#' matrix row names.
#' @param gsets list, list of gene sets. Each gene set is character vector of 
#' Enterz IDs. 
#' The names of the list should match the GeneSet column in the \code{anno} 
#' argument (if it is provided).
#' @param anno list, dataframe that stores a detailed annotation for each  gene 
#' set. 
#' Some of its fields can be ID, GeneSet, PubMed, URLs, etc. The GeneSet field 
#' is mandatory and
#' should have the same names as the \code{gsets}' names. 
#' @param label character,a unique id that identifies the collection of gene 
#' sets
#' @param name character,the collection name to be used in the EGSEA report
#' @param species character, determine the organism of selected gene sets: 
#' "human", "mouse" or "rat".
#' @param min.size integer, the minium number of genes required in a testing 
#' gene set 
#'
#' @return indexed gene set annotation that can be used with other functions in 
#' the package.
#' Each annotation is a list of seven elements: \code{original} stores the 
#' original gene sets, 
#' \code{idx} stores the indexed gene sets,  \code{anno} that stores detailed 
#' annotation for each 
#' gene set, \code{label} a unique id that identifies the collection of gene 
#' sets, 
#' \code{featureIDs} stores the entrezIDs used in building the annotation, 
#' \code{species}
#'  stores that organism name of gene sets and \code{name}  stores the 
#' collection name 
#' to be used in the EGSEA report.
#' 
#' @importFrom limma ids2indices
#' @import EGSEAdata
#' @export 
#' @examples
#' library(EGSEAdata) 
#' data(il13.data)
#' v = il13.data$voom
#' kegg = buildIdxEZID(entrezIDs=rownames(v$E), species="human", 
#' msigdb.gsets="none", 
#'          kegg.updated=FALSE, kegg.exclude = c("Metabolism"))
#' gsets = kegg$kegg$original[1:50]
#' gs.annots = buildCustomIdxEZID(entrezIDs=rownames(v$E), gsets= gsets, 
#' species="human")
#' names(gs.annots)
#' 
#' 

buildCustomIdxEZID <- function(entrezIDs, gsets, anno=NULL,label="custom", 
name="Custom",
        species="Human", min.size=1){   
    entrezIDs = as.character(entrezIDs)

    print("Building custom pathways annotation object ... ")
    for (j in 1:length(gsets)){
        gsets[[j]] = as.character(gsets[[j]])
    }       
    gsets.ez = gsets
    gsets.idx = ids2indices(gsets.ez, entrezIDs, remove.empty=TRUE)
    gsets.ez = gsets.ez[names(gsets.idx)]            
    gsets.names = names(gsets.ez)
    names(gsets.idx) = gsets.names
    names(gsets.ez) = gsets.names   
    
    if (is.null(anno)){
        anno = data.frame(ID=paste0(label, seq(1, length(gsets.idx))), 
GeneSet=gsets.names)
        rownames(anno) = gsets.names
    }else{
        anno = anno[match(gsets.names, anno[, "GeneSet"]), ]        
    }
    anno[, "NumGenes"] = paste0(sapply(gsets.idx, length), 
'/',sapply(gsets.ez, length))
    if (label %in% c("h", "c1", "c2", "c3", "c4", "c5", "c6", "c7")){   
        
        name = msigdb.gs.names[[label]]
    }
        
    species = tolower(species)
    if (species %in% human.names){
        species = "Homo sapiens"
    }else if (species %in% mouse.names){
        species = "Mus musculus"
    }
    else if (species %in% rat.names){
        species = "Rattus norvegicus"
    }
    
    gs.annot = list(original=gsets.ez, idx=gsets.idx, anno=anno, 
label=label, featureIDs=entrezIDs,
            species=species, name=name)
    gs.annot = getGsetAnnot(gs.annot=gs.annot, min.size=min.size)   
    if (length(gs.annot$idx) == 0)
        print(paste0("The cutsom gene set collection is empty."))
    return(gs.annot)
}



#' @title Generate Gene Set Collection Indexes from the MSigDB and KEGG 
#' Databases
#' 
#' @description It prepares the MSigDB and KEGG gene set collections to be used 
#' for the EGSEA analysis. 
#' 
#' @details It indexes the MSigDB and KEGG gene sets and loads gene set 
#' annotation.
#'   
#' @param entrezIDs character, a vector that stores the Entrez Gene IDs tagged 
#' in your dataset. 
#' The order of the Entrez Gene IDs should match those of the count/expression 
#' matrix row names.
#' @param species character, determine the organism of selected gene sets: 
#' "human", "mouse" or "rat".
#' @param msigdb.gsets character, a vector determines which gene set 
#' collections should be used.
#' It can take values from this list: "h", "c1", "c2", "c3", "c4", "c5", 
#' "c6","c7". "h" and "c1"
#' are human specific. If NULL, all available gene set collections are loaded. 
#' If "none", 
#' MSigDB collections are excluded. 
#' @param kegg.updated logical, set to TRUE if you want to download the most 
#' recent KEGG pathways.
#' @param kegg.exclude character, vector used to exclude KEGG pathways of 
#' specific type(s): 
#' Disease, Metabolism, Signaling. If "all", none fo the KEGG collections is 
#' included.  
#' @param min.size integer, the minium number of genes required in a testing 
#' gene set
#' @param rdata.dir character, directory from which the MSigDB collections are 
#' loaded. 
#' If NULL, EGSEA tries to load the data from \pkg{EGSEAdata}. 
#'
#' @return indexed gene set annotation that can be used with other functions in 
#' the package.
#' Each annotation is a list of seven elements: \code{original} stores the 
#' original gene sets, 
#' \code{idx} stores the indexed gene sets,  \code{anno} that stores detailed 
#' annotation for each 
#' gene set, \code{label} a unique id that identifies the collection of gene 
#' sets, 
#' \code{featureIDs} stores the entrezIDs used in building the annotation, 
#' \code{species}
#'  stores that organism name of gene sets and \code{name}  stores the 
#' collection name 
#' to be used in the EGSEA report.
#' 
#' @import EGSEAdata
#' @export 
#' @examples 
#' library(EGSEAdata)
#' data(il13.data)
#' v = il13.data$voom
#' gs.annots = buildIdxEZID(entrezIDs=rownames(v$E), species="human",
#'          msigdb.gsets = c("h", "c2"),
#'          kegg.exclude = c("Metabolism"))
#' names(gs.annots)



buildIdxEZID <- function(entrezIDs, species="human", 
        msigdb.gsets=NULL, kegg.updated=FALSE, kegg.exclude=c(), 
min.size = 1,
        rdata.dir=NULL ){
    if (length(msigdb.gsets) == 1 && tolower(msigdb.gsets[1]) == "none")
        gs.annots = list()
    else
        gs.annots = buildMSigDBIdxEZID(entrezIDs=entrezIDs, 
            geneSets=msigdb.gsets,
            species = species,  min.size = min.size,
            rdata.dir=rdata.dir) # type ?buildMSigDBIdxEZID to see 
#docs! 
    if (length(kegg.exclude) == 1 && tolower(kegg.exclude[1]) == "all")
        return(gs.annots)
    gs.annot = buildKEGGIdxEZID(entrezIDs=entrezIDs,species = species,  
min.size= min.size,
            rdata.dir=rdata.dir, updated = kegg.updated) # type 
#?buildKEGGIdxEZID to see docs! 
    ### shall you want to exclude the KEGG metabolic pathways
    sel = ! tolower(gs.annot$anno[, "Type"]) %in% tolower(kegg.exclude)
    gs.annot$idx = gs.annot$idx[sel]
    gs.annot$original = gs.annot$original[sel]
    gs.annot$anno = gs.annot$anno[sel, ]
    ### Otherwise, continue here
    gs.annots[["kegg"]] = gs.annot
    rm(gs.annot)
    return(gs.annots)
}

