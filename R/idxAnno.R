#Ensemble of Gene Set Enrichment Analyses 
# 
# Author: Monther Alhamdoosh, E:m.hamdoosh@gmail.com
###############################################################################


species.fullToShort = list()
species.fullToShort[["homo sapiens"]] = "human"
species.fullToShort[["mus musculus"]] = "mouse"
species.fullToShort[["rattus norvegicus"]] = "rat"
# "Homo sapiens", "Mus musculus", "Rattus norvegicus",
# "Danio rerio", "Macaca mulatta"

loadGeneSetDBCategoryLabels <- function(){
    genesetdb.gs.labels = list("GeneSetDB Drug/Chemical"="gsdbdrug",   
            "GeneSetDB Disease/Phenotype" = "gsdbdis",
            "GeneSetDB Gene Ontology" = "gsdbgo",
            "GeneSetDB Pathway" = "gsdbpath",
            "GeneSetDB Gene Regulation" = "gsdbreg")
    return(genesetdb.gs.labels)
}

normalizeSpecies <- function(species){
    human.names = c("human", "homo sapiens", "hs")
    mouse.names = c("mouse", "mus musculus", "mm")
    rat.names = c("rat", "rattus norvegicus" , "rn")
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
    return(species)
}

loadKeggData <- function(){
    data("kegg.pathways", package="EGSEAdata")
    if (is.null(kegg.pathways))
        stop("Failed to load the KEGG pathway data.")
    return (kegg.pathways)
}


#' @title Functions to create gene set collection indexes for EGSEA
#' 
#' @description \code{buildIdx} indexes the MSigDB, KEGG  
#' and GeneSetDB collections to be used for the EGSEA analysis. 
#' 
#' @details \code{buildIdx} indexes the MSigDB, KEGG 
#' and GeneSetDB gene set collections, and loads gene set annotation.
#'   
#' @param entrezIDs character, a vector that stores the Entrez Gene IDs tagged 
#' in your dataset. 
#' The order of the Entrez Gene IDs should match those of the count/expression 
#' matrix row names.
#' @param species character, determine the organism of selected gene sets: 
#' "human", "mouse" or "rat".
#' @param msigdb.gsets character, a vector determines which gene set 
#' collections should be used from MSigDB.
#' It can take values from this list: "h", "c1", "c2", "c3", "c4", "c5", 
#' "c6","c7". "h" and "c1"
#' are human specific. If "all", all available gene set collections are loaded. 
#' If "none", 
#' MSigDB collections are excluded.
#' @param gsdb.gsets character, a vector determines which gene set collections
#' are loaded from the GeneSetDB.
#' It takes "none", "all", "gsdbdis", "gsdbgo", "gsdbdrug", "gsdbpath" or "gsdbreg". 
#' "none" excludes the GeneSetDB collections.
#' "all" includes all the GeneSetDB collections.  
#' "gsdbdis" to load the disease collection, "gsdbgo" to load the GO terms collection,
#'  "gsdbdrug" to load the drug/chemical collection, 
#' "gsdbpath" to load the pathways collection and "gsdbreg" to load the gene regulation
#' collection. 
#' @param go.part logical, whether to partition the GO term collections into the
#' three GO domains: CC, MF and BP or use the entire collection all together. 
#' @param kegg.updated logical, set to TRUE if you want to download the most 
#' recent KEGG pathways.
#' @param kegg.exclude character, vector used to exclude KEGG pathways of 
#' specific type(s): 
#' Disease, Metabolism, Signaling. If "all", none fo the KEGG collections is 
#' included.  
#' @param min.size integer, the minium number of genes required in a testing 
#' gene set
#'
#' @return \code{buildIdx} returns a list of gene set collection indexes, where
#' each element of the list is an object of the class GSCollectionIndex. 
#' 
#' @import EGSEAdata
#' @export 
#' 
#' @name buildIdx
#' @aliases buildIdx,egsea-index
#' @rdname egsea-index
#' 
#' @examples 
#' # example of buildIdx
#' library(EGSEAdata)
#' data(il13.data)
#' v = il13.data$voom
#' gs.annots = buildIdx(entrezIDs=rownames(v$E), species="human",
#'          msigdb.gsets = c("h", "c5"),
#' 			go.part = TRUE,
#'          kegg.exclude = c("Metabolism"))
#' names(gs.annots)

buildIdx <- function(entrezIDs, species="human", 
        msigdb.gsets="all",        
        gsdb.gsets = "none",
        go.part = FALSE,
        kegg.updated=FALSE, kegg.exclude=c(), 
        min.size = 1){
    if (length(msigdb.gsets) == 1 && tolower(msigdb.gsets[1]) == "none")
        gs.annots = list()
    else
        gs.annots = buildMSigDBIdx(entrezIDs=entrezIDs,            
                species = species,  
                geneSets = msigdb.gsets,
                go.part = go.part,
                min.size = min.size)
    if (!is.null(gsdb.gsets) && tolower(gsdb.gsets[1]) != "none")
        gs.annots = c(gs.annots, buildGeneSetDBIdx(entrezIDs = entrezIDs, 
                        species = species, geneSets = gsdb.gsets, 
                        go.part = go.part,
                        min.size=min.size))
    if (length(kegg.exclude) == 1 && tolower(kegg.exclude[1]) == "all")
        return(gs.annots)
    gs.annot = buildKEGGIdx(entrezIDs=entrezIDs,species = species,  
            min.size= min.size, updated = kegg.updated, exclude=kegg.exclude) 
    gs.annots[["kegg"]] = gs.annot
    rm(gs.annot)
    return(gs.annots)
}


#' @title Gene Set Collection Index from the KEGG Database
#' 
#' @description \code{buildKEGGIdx} prepares the KEGG pathway collection to 
#' be used for the EGSEA analysis. 
#' 
#' @details \code{buildKEGGIdx} indexes the KEGG pathway gene sets and 
#' loads gene set annotation.
#'   
#' @inheritParams entrezIDs 
#' @inheritParams species 
#' @inheritParams min.size 
#' @param updated logical, set to TRUE if you want to download the most recent 
#' KEGG pathways.
#' @param exclude character, vector used to exclude KEGG pathways of 
#' specific category. Accepted values are "Disease", "Metabolism", or "Signaling".
#'
#' @return \code{buildKEGGIdx} returns an object of the class GSCollectionIndex. 
#' 
#' @import EGSEAdata
#' @export 
#' 
#' @name buildKEGGIdx
#' @aliases buildKEGGIdx,egsea-index
#' @rdname egsea-index
#' 
#' @examples 
#' # example of buildKEGGIdx
#' library(EGSEAdata)
#' data(il13.data)
#' v = il13.data$voom
#' gs.annots = buildKEGGIdx(entrezIDs=rownames(v$E), species="human")
#' 

buildKEGGIdx <- function(entrezIDs, species = "human", min.size=1, 
          updated=FALSE, exclude = c()) {   
    species = normalizeSpecies(species)
    updatedSuccess = FALSE
    entrezIDs = as.character(entrezIDs)    
    kegg = NULL
    print("Building KEGG pathways annotation object ... ")      
    if (!updated){
        kegg.pathways = loadKeggData()
        kegg = kegg.pathways[[species.fullToShort[[tolower(species)]]]]
    }else{
        tryCatch({ 
            kegg <- kegg.gsets(species = 
            species.fullToShort[[tolower(species)]], id.type = "kegg")
            updatedSuccess <- TRUE
        },
        error = function(e){    
            warning("KEGG pathways have not been updated successfully.")                
            kegg.pathways <- loadKeggData() 
            kegg <- kegg.pathways[[species.fullToShort[[tolower(species)]]]]
        })
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
    if (!updatedSuccess){
        db.info = egsea.data(simple=TRUE, returnInfo=TRUE)
        ver = db.info$kegg$info$version
        dat = db.info$kegg$info$date
    }else{
        ver = "NA"
        dat = date()
    }
    
    gs.annot = GSCollectionIndex(original = gsets.ez,
            idx = gsets,
            anno = anno,                
            featureIDs = entrezIDs,
            species = species,
            name = "KEGG Pathways",
            label = "kegg",
            version = ver,
            date = dat)
    
    gs.annot = selectGeneSets(gs.annot, min.size=min.size)   
    if (length(gs.annot@idx) == 0)
        cat("KEGG pathway collection is empty.\n")
    ### shall you want to exclude the KEGG metabolic pathways
    if (length(exclude) > 0){
        sel = ! tolower(gs.annot@anno[, "Type"]) %in% tolower(exclude)
        gs.annot@idx = gs.annot@idx[sel]
        gs.annot@original = gs.annot@original[sel]
        gs.annot@anno = gs.annot@anno[sel, ]
    }
   
    return(gs.annot)
}


#' @title Gene Set Collection Indexes from the MSigDB Database
#' 
#' @description \code{buildMSigDBIdx} prepares the MSigDB gene set collections
#'  to be used for the EGSEA analysis. 
#' 
#' @details \code{buildMSigDBIdx} indexes the MSigDB gene sets and loads gene 
#' set annotation.
#'   
#' @inheritParams entrezIDs 
#' @inheritParams species 
#' @param geneSets character, a vector determines which gene set collections 
#' should be used. For MSigDB, it can take values from this list: 
#' "all", "h", "c1", "c2", "c3", "c4", "c5", "c6","c7". "c1"
#' is human specific. For GeneSetDB, it takes 
#' "all", "gsdbdis", "gsdbgo", "gsdbdrug", "gsdbpath" or "gsdbreg".
#' "gsdbdis" is to load the disease collection, "gsdbgo" to load the GO terms collection,
#'  "gsdbdrug" to load the drug/chemical collection, 
#' "gsdbpath" to load the pathways collection and "gsdbreg" to load the gene regulation
#' collection.  If "all", all available gene set collections are loaded. 
#' @inheritParams go.part 
#' @inheritParams min.size 
#'
#' @return \code{buildMSigDBIdx} returns a list of gene set collection indexes, where
#' each element of the list is an object of the class GSCollectionIndex. 
#' 
#' @import EGSEAdata
#' @export 
#' 
#' @name buildMSigDBIdx
#' @aliases buildMSigDBIdx,egsea-index
#' @rdname egsea-index
#' 
#' @examples 
#' # example of buildMSigDBIdx
#' library(EGSEAdata)
#' data(il13.data)
#' v = il13.data$voom
#' gs.annots = buildMSigDBIdx(entrezIDs=rownames(v$E), species="human",
#' geneSets=c("h", "c2"))
#' names(gs.annots)


buildMSigDBIdx <- function(entrezIDs, 
        species="Homo sapiens", geneSets="all", 
        go.part = FALSE, min.size=1){
    geneSets = tolower(geneSets)
    stopifnot(geneSets %in% c("all", "h", "c1", "c2", "c3", "c4", "c5", 
                "c6","c7"))
    msigdb.gs.names = list(c1="c1 Positional Gene Sets", 
        c2="c2 Curated Gene Sets", 
        c3="c3 Motif Gene Sets", 
        c4="c4 Computational Gene Sets", c5="c5 GO Gene Sets", 
        c6="c6 Oncogenic Signatures", c7="c7 Immunologic Signatures",
        h="h Hallmark Signatures")
    #stopifnot(!is.null(geneSets))   
    species = normalizeSpecies(species)
    if (species == "Homo sapiens"){        
        if (length(geneSets) == 1 && geneSets == "all")
            geneSets = c("h", "c1", "c2", "c3", "c4", "c5", 
                "c6","c7")
    }else if (species == "Mus musculus"){       
        if (length(geneSets) == 1 && geneSets == "all")
            geneSets = c("h", "c2", "c3", "c4", "c5", "c6","c7")
    }
    else{
        stop("Unrecognized species. MSigDB supports human and mouse 
			only.")
    }
    
    entrezIDs = as.character(entrezIDs)
    gs.annots = list() # Gene set annotation indexes
    # Read all gene set collections
    print("Loading MSigDB Gene Sets ... ")
    # SymbolIdentifier(), EntrezIdentifier()   
    data("msigdb", package="EGSEAdata")
    gsc.all = msigdb              
    if (is.null(gsc.all))
        stop("Failed to load the MSigDB gene set collection data")
    gsc.all = gsc.all[names(gsc.all) == "GENESET"]          
    names(gsc.all) = sapply(gsc.all, function(x) x["STANDARD_NAME"])
    organisms = sapply(gsc.all, function(x) x["ORGANISM"])
    if (species == "Homo sapiens")
        gsc.all = gsc.all[organisms == "Homo sapiens"] #  "Homo 
#sapiens","Mus musculus", "Rattus norvegicus","Danio rerio","Macaca mulatta"
    types = tolower(sapply(gsc.all, function(x) x["CATEGORY_CODE"])) 
    db.info = egsea.data(simple=TRUE, returnInfo=TRUE)
#character vector c1 c2 c3 c4 or c5
    ## process each gene set collection
    for (geneSet in geneSets){     
        gs.annot = GSCollectionIndex(version = db.info$msigdb$info$version,
                date = db.info$msigdb$info$date)        
        gsc.small = gsc.all[types == geneSet] 
        if (species == "Mus musculus"){
            if (geneSet %in% c("h", "c2", "c3", "c4", "c5", 
                        "c6","c7")){                
                geneSet1 = ifelse(geneSet == "h", "H", geneSet)               
                data(list=paste0("Mm.", geneSet1), package="EGSEAdata")                
                gs.annot@original =get(paste0("Mm.", geneSet1))
            }
            else{
                warning(paste0("Unsupported gene set collection for Mus", 
					" Musculus ... ", geneSet))
                next
            }
        }
        else if (species == "Homo sapiens")
            gs.annot@original = sapply(gsc.small, function (x) 
                strsplit(x["MEMBERS_EZID"], ",")[[1]])      
        print(paste0("Loaded gene sets for the collection ", geneSet,  
					" ..."))
        gs.annot@idx = ids2indices(gs.annot@original, entrezIDs, 
                remove.empty=TRUE) # pathways
        gsc.small = gsc.small[names(gs.annot@idx)] # remove 
        gs.annot@original = 
                gs.annot@original[names(gs.annot@idx)] 
#empty sets
        print(paste0("Indexed the collection ", geneSet, 
                        " ..."))  
        if (length(gs.annot@idx) == 0){
            cat(paste0("None of the genes in ", geneSet, 
                            " are mapped to your gene IDs\n"))          
        }        
        # create annotation frame
        if (length(gs.annot@idx) > 0){                       
            tmp = rep(NA, length(gs.annot@idx)) 
            gs.annot@anno = data.frame(ID=tmp, 
                GeneSet=names(gs.annot@idx), BroadUrl=tmp, Description=tmp,
                PubMedID=tmp, NumGenes=rep(0, 
                length(gs.annot@idx)), Contributor=tmp)
            gs.annot@anno[, "ID"] = sapply(gsc.small, function(x) 
                x["SYSTEMATIC_NAME"])
            gs.annot@anno[, "BroadUrl"] = 
                paste0("http://www.broadinstitute.org/gsea/msigdb/cards/", 
                names(gs.annot@idx), 
                ".html")    
            gs.annot@anno[, "Description"] = sapply(gsc.small, 
                function(x) x["DESCRIPTION_BRIEF"])
            gs.annot@anno[, "PubMedID"] = sapply(gsc.small, 
                function(x) x["PMID"])
            gs.annot@anno[, "NumGenes"] = 
                paste0(sapply(gs.annot@idx, length), '/',
                        sapply(gs.annot@original, length)) 
    
            gs.annot@anno[, "Contributor"] = sapply(gsc.small, 
                function(x) x["CONTRIBUTOR"])   # KEGG      
            rownames(gs.annot@anno) = names(gs.annot@original)
            print(paste0("Created annotation for the collection ", geneSet, 
                " ..."))
            gs.annot = selectGeneSets(gs.annot, min.size=min.size)   
            if (length(gs.annot@idx) == 0)
                cat(paste0("MSigDB ", gs.annot@label, " gene set ", 
                                "collection is empty.\n"))
        }        
        gs.annot@label = geneSet    
        gs.annot@featureIDs = entrezIDs
        gs.annot@species = species
        gs.annot@name = msigdb.gs.names[[gs.annot@label]]     
        
        if (geneSet == "c5" && length(gs.annot@idx) > 0){
            external.urls = sapply(gsc.small, 
                    function(x) x["EXTERNAL_DETAILS_URL"])
            goID = gsub(".*(GO:[0-9]+).*$","\\1", 
                    external.urls)
            ontology = sapply(gsc.small, 
                    function(x) x["SUB_CATEGORY_CODE"])
            gs.annot@anno = cbind(gs.annot@anno, Ontology=ontology, 
                GOID=goID)
            gs.annot@anno = droplevels(gs.annot@anno)
            if (go.part){
                new.labels = c()
                for (domain in c("BP", "CC", "MF")){                
                    gs.annot.tmp = selectGeneSets(gs.annot,
                            gs.names = gs.annot@anno[
                                    gs.annot@anno[, "Ontology"] == domain ,
                                    "GeneSet"])
                    gs.annot.tmp@label = paste0(gs.annot@label, domain)
                    gs.annot.tmp@name = paste0(gs.annot@name, " (", domain, ")")
                    gs.annots[[gs.annot.tmp@label]] = gs.annot.tmp
                    if (length(gs.annot.tmp@idx) == 0)
                        cat(paste0("MSigDB ", gs.annot.tmp@label, " gene set ", 
                                        "collection is empty.\n"))  
                    new.labels = c(new.labels, gs.annot.tmp@label)
                }
                cat(paste0("MSigDB ", gs.annot@label, " gene set ", 
                                "collection has been partitioned into \n",
                                paste(new.labels, collapse=", "), "\n"))
            }else
                gs.annots[[geneSet]] = gs.annot
        }else
            gs.annots[[geneSet]] = gs.annot
    }
    return(gs.annots)
    
}


#' @title Gene Set Collection Indexes from the GeneSetDB Database
#' 
#' @description \code{buildGeneSetDBIdx} prepares the GeneSetDB gene set 
#' collections to be used for the EGSEA analysis. 
#' 
#' @details \code{buildGeneSetDBIdx} indexes the GeneSetDB gene sets and 
#' loads gene set annotation.
#'   
#' @inheritParams entrezIDs 
#' @inheritParams species 
#' @inheritParams geneSets
#' @inheritParams go.part 
#' @inheritParams min.size 
#'
#' @return \code{buildGeneSetDBIdx} returns a list of gene set collection indexes, where
#' each element of the list is an object of the class GSCollectionIndex. 
#' 
#' @import EGSEAdata
#' @export 
#' 
#' @name buildGeneSetDBIdx
#' @aliases buildGeneSetDBIdx,egsea-index
#' @rdname egsea-index
#' 
#' @examples 
#' # example of buildGeneSetDBIdx
#' library(EGSEAdata)
#' data(il13.data)
#' v = il13.data$voom
#' gs.annots = buildGeneSetDBIdx(entrezIDs=rownames(v$E), species="human")
#' names(gs.annots)
#' 
#' 


buildGeneSetDBIdx <- function(entrezIDs, species, geneSets="all", 
        go.part = FALSE, min.size=1){
    geneSets = tolower(geneSets)
    stopifnot(geneSets %in% c("all", "gsdbdis", "gsdbgo", "gsdbdrug", 
                    "gsdbpath" , "gsdbreg"))
    print("Loading GeneSetDB Gene Sets ... ")
    species = normalizeSpecies(species)     
    data(list=paste0('gsetdb.', 
            species.fullToShort[[tolower(species)]]), 
            package="EGSEAdata")
    gsetdb.all = get(paste0('gsetdb.', 
                species.fullToShort[[tolower(species)]]))
    if (is.null(gsetdb.all)){
        stop("Failed to load the GeneSetDB collection data")
    }
    
    gs.annot = buildCustomIdx(geneIDs = entrezIDs, gsets = 
        gsetdb.all$original,
        anno = gsetdb.all$anno, label = "gsDB", 
        name="GeneSetDB Gene Sets", species = species)
    db.info = egsea.data(simple=TRUE, returnInfo=TRUE)
    gs.annot@version = db.info$gsetdb$info$version
    gs.annot@date = db.info$gsetdb$info$date
    
    categories = unique(as.character(gs.annot@anno$Category))
    gs.annots = list()
    genesetdb.gs.labels = loadGeneSetDBCategoryLabels()
    for (cat in categories){
        label = genesetdb.gs.labels[[cat]]
        gs.annot.cat = selectGeneSets(gs.annot, 
            gs.names = as.character(gs.annot@anno[
                gs.annot@anno[,"Category"] == cat,"GeneSet"]),  
            min.size=min.size)
        if (length(gs.annot.cat@idx) == 0){
            cat(paste0("GeneSetDB ", label," gene set 
				collection is empty.\n"))             
        }
       
        gs.annot.cat@label = label
        gs.annot.cat@name = cat
        gs.annot.cat@anno = gs.annot.cat@anno[, 
            -c(match("Category", names(gs.annot@anno)))]
        if (label == "gsdbgo" && length(gs.annot.cat@idx) > 0){
            go.terms = gsub(".*(GO:[0-9]+).*$","\\1", 
                    gs.annot.cat@anno[, "GeneSet"]) 
            xx = GOTERM                
            sel.gsets = character(0)
            goID = character(0)
            ontology = character(0)
            no.term = 0 # no GO term found
            for (i in 1:length(go.terms)){              
                temp = xx[[go.terms[i]]]
                if (!is.null(temp)){                        
                    sel.gsets = c (sel.gsets, 
                            as.character(gs.annot.cat@anno[i, "GeneSet"]))
                    ontology = c(ontology, Ontology(temp))
                    goID = c(goID, go.terms[i])
                }else{
                    no.term = no.term + 1
                }
            }            
            if (no.term > 0){
                cat(paste0(no.term, " gene sets from the GeneSetDB ",
                    label, " collection do not have valid GO ID.\n",
                    "They will be removed. \n"))
                gs.annot.cat = selectGeneSets(gs.annot.cat, sel.gsets)
            }            
            gs.annot.cat@anno = cbind(gs.annot.cat@anno,  
                    Ontology = ontology, GOID = goID)
            gs.annot.cat@anno = droplevels(gs.annot.cat@anno) 
#            print(head(gs.annot.cat@anno))
            if (go.part){
                new.labels = c()
                for (domain in c("BP", "CC", "MF")){                
                    gs.annot.tmp = selectGeneSets(gs.annot.cat,
                            gs.names = gs.annot.cat@anno[
                                    gs.annot.cat@anno[, "Ontology"] == domain ,
                                    "GeneSet"])
                    gs.annot.tmp@label = paste0(gs.annot.cat@label, domain)
                    gs.annot.tmp@name = paste0(gs.annot.cat@name, " (", domain, ")")
                    gs.annots[[gs.annot.tmp@label]] = gs.annot.tmp
                    if (length(gs.annot.tmp@idx) == 0)
                        cat(paste0("GeneSetDB ", gs.annot.tmp@label, " gene set ", 
                                        "collection is empty.\n"))
                    new.labels = c(new.labels, gs.annot.tmp@label)
                }
                cat(paste0("GeneSetDB ", gs.annot.cat@label, " gene set ", 
                                "collection has been partitioned into \n",
                                paste(new.labels, collapse=", "), "\n"))
            }else
                gs.annots[[label]] = gs.annot.cat
        }else
            gs.annots[[label]] = gs.annot.cat
    }
    if (length(geneSets) == 1 && geneSets == "all")        
        return(gs.annots)
    else
        return(gs.annots[geneSets])
}




#' @title Custom Gene Set Collection Index 
#' 
#' @description \code{buildCustomIdx} creates a gene set collection from a 
#' given list of gene sets to be used for the EGSEA analysis. 
#' 
#' @details \code{buildCustomIdx} indexes newly created gene sets and 
#' attach gene set annotation if provided.
#'   
#' @param geneIDs character, a vector that stores the Gene IDs 
#' tagged in your dataset. The order of the Gene IDs must match those of 
#' the count/expression matrix row names. Gene IDs can be in any annotation, e.g.,
#' Symbols, Ensembl, etc., as soon as the parameter \code{gsets} uses the same
#' Gene ID annotation.  
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
#' @inheritParams species 
#' @inheritParams min.size 
#'
#' @return \code{buildCustomIdx} returns an object of the class GSCollectionIndex.
#' 
#' @importFrom limma ids2indices
#' @import EGSEAdata
#' @export 
#' 
#' @name buildCustomIdx
#' @aliases buildCustomIdx,egsea-index
#' @rdname egsea-index
#' 
#' @examples
#' # example of buildCustomIdx
#' library(EGSEAdata) 
#' data(il13.data)
#' v = il13.data$voom
#' data(kegg.pathways)
#' gsets = kegg.pathways$human$kg.sets[1:50]
#' gs.annot = buildCustomIdx(geneIDs=rownames(v$E), gsets= gsets, 
#' species="human")
#' class(gs.annot)
#' 
#' 

buildCustomIdx <- function(geneIDs, gsets, anno=NULL,label="custom", 
        name="User-Defined Gene Sets", species="Human", min.size=1){   
    geneIDs = as.character(geneIDs)
    for (j in 1:length(gsets)){
        gsets[[j]] = as.character(gsets[[j]])
    }       
    gsets.ez = gsets
    gsets.idx = ids2indices(gsets.ez, geneIDs, remove.empty=TRUE)
    gsets.ez = gsets.ez[names(gsets.idx)]            
    gsets.names = names(gsets.ez)
    names(gsets.idx) = gsets.names
    names(gsets.ez) = gsets.names   
    
    if (is.null(anno)){
        anno = data.frame(ID=paste0(label, seq(1, length(gsets.idx))), 
            GeneSet=gsets.names)
    }else{
        stopifnot("GeneSet" %in% colnames(anno))
        anno = anno[match(gsets.names, anno[, "GeneSet"]), ]        
    }
    rownames(anno) = gsets.names
    anno[, "NumGenes"] = paste0(sapply(gsets.idx, length), 
        '/',sapply(gsets.ez, length))

    try(species <- normalizeSpecies(species), silent = TRUE)
    gs.annot = GSCollectionIndex(original = gsets.ez,
            idx = gsets.idx,
            anno = anno,                
            featureIDs = geneIDs,
            species =species,
            name = name,
            label = label,
            version = "NA",
            date = date())
    gs.annot = selectGeneSets(gs.annot, min.size=min.size)   
    if (length(gs.annot@idx) == 0)
        cat(paste0("The cutsom gene set collection is empty.\n"))
    print(paste0("Created the ", name, " collection ... "))
    return(gs.annot)
}

#' @title Gene Set Collection Index from a GMT file
#' 
#' @description \code{buildGMTIdx} creates a gene set collection from a 
#' given GMT file to be used for the EGSEA analysis. 
#' 
#' @details \code{buildGMTIdx} indexes newly created gene sets and 
#' attach gene set annotation if provided.
#'   
#' @inheritParams geneIDs
#' @param gmt.file character, the path and name of the GMT file
#' @param anno.cols integer, number of columns in the GMT file that are 
#' used for annotation. These columns should be inserted immediately after
#' the second column. 
#' @param anno.col.names character, vector of the names of the annotation
#' columns. 
#' @inheritParams label 
#' @inheritParams name 
#' @inheritParams species 
#' @inheritParams min.size 
#'
#' @return \code{buildGMTIdx} returns an object of the class GSCollectionIndex.
#' 
#' @export 
#' 
#' @name buildGMTIdx
#' @aliases buildGMTIdx,egsea-index
#' @rdname egsea-index
#' 
#' @examples
#' # example of buildGMTIdx
#' library(EGSEAdata) 
#' data(il13.data)
#' v = il13.data$voom
#' #gs.annot = buildGMTIdx(geneIDs=rownames(v$E), gsets= gmt.file, 
#' #species="human")
#' #class(gs.annot)
#' 
#' 
buildGMTIdx <- function(geneIDs, gmt.file, anno.cols = 0, 
        anno.col.names = NULL, label="gmtcustom", 
        name="User-Defined GMT Gene Sets", species="Human", min.size=1){
    fc <- file(gmt.file)
    gsets.raw <- strsplit(readLines(fc), "\t")
    close(fc)
    
    gsets = lapply(gsets.raw, function(x) x[-(1:(anno.cols + 2))])
    gsets.names = sapply(gsets.raw, function(x) x[1])
    desc = sapply(gsets.raw, function(x) x[2])
    names(gsets) = gsets.names
    anno = data.frame(ID=paste0(label, seq(1, length(gsets.names))), 
            GeneSet=gsets.names, Description=desc)
    if (anno.cols > 0){
        if (is.null(anno.col.names))
            anno.col.names = paste0("Anno", seq(1, anno.cols))        
        stopifnot(length(anno.col.names) == anno.cols)
        
        for (cl in 3:(3+anno.cols-1)){
            annocol = sapply(gsets.raw, function(x) x[cl])
            anno = cbind(anno, annocol)        
        }        
        colnames(anno)[4:(4+anno.cols-1)] = anno.col.names        
    }
    rownames(anno) = gsets.names
    
    return(buildCustomIdx(geneIDs, gsets,
                    anno, label, name, species, min.size))
}
