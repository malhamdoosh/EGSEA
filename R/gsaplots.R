#Ensemble of Gene Set Enrichment Analyses 
# 
# Author: Monther Alhamdoosh, E:m.hamdoosh@gmail.com
###############################################################################
# Plot GO graphs using topGO package

plotGOGraphs <- function(egsea.results, gene.ezids, gs.annot, gsa.dir, sort.by){
    print("GO graphs are being generated for top-ranked GO terms based on 
p-values ... ")

    
    gl = rep(0, length(gs.annot$featureIDs))
    names(gl) = gs.annot$featureIDs
    go.dir = paste0(gsa.dir, "/go-graphs/") 
    dir.create(file.path(go.dir), showWarnings = FALSE)
    contrast.names = names(egsea.results)
    file.name = paste0(go.dir, sub(" - ", "-", contrast.names), "-", 
gs.annot$label, "-top-")
    if (file.exists(paste0(file.name[length(file.name)], "CC.png")))
        return()
    if (tolower(gs.annot$species) %in% c("human", "homo sapiens")){
        mappingDB = "org.Hs.eg.db"
    }else if (tolower(gs.annot$species) %in% c("mouse", "mus musculus")){
        mappingDB = "org.Mm.eg.db"
    }
#   print(mappingDB)
    capture.output(topGOdataBP <- new("topGOdata",ontology = "BP", allGenes = gl,
            geneSel = topDiffGenes, nodeSize = 10, annot = 
annFUN.org, mapping=mappingDB))
capture.output(topGOdataMF <- new("topGOdata",ontology = "MF", allGenes = gl,
            geneSel = topDiffGenes, nodeSize = 10, annot = 
annFUN.org, mapping=mappingDB))
capture.output(topGOdataCC <- new("topGOdata",ontology = "CC", allGenes = gl,
            geneSel = topDiffGenes, nodeSize = 10, annot = 
annFUN.org, mapping=mappingDB))
#   print(topGOdataBP)
#   print(summary(topGOdataBP))
    
    go.subsets = list() 
    go.ids = as.character(gs.annot$anno[match(rownames(egsea.results[[1]]), 
                    gs.annot$anno[, "GeneSet"]), "GOID"])

    go.subsets[["BP"]] = go.ids[go.ids %in% topGOdataBP@graph@nodes]
    go.subsets[["MF"]] = go.ids[go.ids %in% topGOdataMF@graph@nodes]
    go.subsets[["CC"]] = go.ids[go.ids %in% topGOdataCC@graph@nodes]    
#   print(go.subsets)
    for (i in 1:length(contrast.names)){            
        print(contrast.names[i])
        scores = egsea.results[[i]][, sort.by]
        max.score = 1.001
        noSig = 5
        if (max(scores) > 1)
            scores = (scores - min(scores)) / (max(scores) - 
min(scores)) 
        scores = scores + 0.001
        names(scores) = 
as.character(gs.annot$anno[match(rownames(egsea.results[[i]]), 
                                gs.annot$anno[, 
"GeneSet"]), "GOID"])       

        tryCatch({
            scores.sub = rep(max.score, 
length(topGOdataBP@graph@nodes))
            names(scores.sub) = topGOdataBP@graph@nodes
            scores.sub[go.subsets[["BP"]]] = 
scores[go.subsets[["BP"]]]
    #       print(head(pvalues.sub))
    #       print(length(pvalues.sub))
#           print(paste0(file.name[i], "BP.pdf"))
            pdf(paste0(file.name[i], "BP.pdf"))     
            showSigOfNodes(topGOdataBP, scores.sub, 
firstSigNodes=noSig, 
                    useInfo='all', sigForAll=FALSE) # or 
# use printGraph to write out plot to file
            dev.off()
            png(paste0(file.name[i], "BP.png"), width=800, 
height=800)
            showSigOfNodes(topGOdataBP, scores.sub, 
firstSigNodes=noSig, 
                    useInfo='all', sigForAll=FALSE) # or 
# use printGraph to write out plot to file
            dev.off()
        }, error=function(err){
            print(paste("MY_ERROR:  ",err))
            dev.off()
             file.remove(paste0(file.name[i], "BP.pdf"))
             file.remove(paste0(file.name[i], "BP.png"))
            }
        )
#       stop("here")
        
        tryCatch({
            scores.sub = rep(max.score, 
length(topGOdataMF@graph@nodes))
            names(scores.sub) = topGOdataMF@graph@nodes
            scores.sub[go.subsets[["MF"]]] = 
scores[go.subsets[["MF"]]]
            pdf(paste0(file.name[i], "MF.pdf"))
            showSigOfNodes(topGOdataMF, scores.sub, 
firstSigNodes=noSig, 
                    sigForAll=FALSE, useInfo='all') # or 
# use printGraph to write out plot to file
            dev.off()
            png(paste0(file.name[i], "MF.png"), width=800, 
height=800)
            showSigOfNodes(topGOdataMF, scores.sub, 
firstSigNodes=noSig, 
                    sigForAll=FALSE, useInfo='all') # or 
# use printGraph to write out plot to file
            dev.off()
        }, error=function(err){
            print(paste("MY_ERROR:  ",err))
            dev.off()
            file.remove(paste0(file.name[i], "MF.pdf"))
            file.remove(paste0(file.name[i], "MF.png"))
        }
        )       
        
        tryCatch({
            scores.sub = rep(max.score, 
length(topGOdataCC@graph@nodes))
            names(scores.sub) = topGOdataCC@graph@nodes
            scores.sub[go.subsets[["CC"]]] = 
scores[go.subsets[["CC"]]]
#                   print(head(sort(pvalues.sub), 20))
#                   print(pvalues.sub["GO:0005581"])
#                   print(pvalues["GO:0005581"])
            pdf(paste0(file.name[i], "CC.pdf"))
            showSigOfNodes(topGOdataCC, scores.sub, 
firstSigNodes=noSig, 
                    sigForAll=FALSE, useInfo='all') # or 
# use printGraph to write out plot to file
            dev.off()
            png(paste0(file.name[i], "CC.png"), width=800, 
height=800)
            showSigOfNodes(topGOdataCC, scores.sub, 
firstSigNodes=noSig, 
                    sigForAll=FALSE, useInfo='all') # or 
# use printGraph to write out plot to file
            dev.off()
        }, error=function(err){
            print(paste("MY_ERROR:  ",err))
            dev.off()
            file.remove(paste0(file.name[i], "CC.pdf"))
            file.remove(paste0(file.name[i], "CC.png"))
        }
        )
        
    }
    
}

topDiffGenes <- function (allScore) {
    return(allScore < 0.01)
}

# Plot summary plots for each contrast and comparisons

generateSumPlots <- function(egsea.results, baseGSEAs, gs.annot, gsa.dir,        
 
        sum.plot.cutoff=1, sum.plot.axis="p.value"){
    print("Summary plots are being generated ... ")
    
    summary.dir = paste0(gsa.dir, "/summary/")
    dir.create(file.path(summary.dir), showWarnings = FALSE)
    contrast.names = names(egsea.results)   
    
    for(i in 1:length(egsea.results)){      
        file.name = paste0(summary.dir, sub(" - ", "-", 
        contrast.names[i]), "-", gs.annot$label, "-summary")        

        if (!file.exists(paste0(file.name, ".dir.png"))){
            plot.data = generatePlotData(egsea.results[[i]], 
                     gs.annot, sum.plot.cutoff, sum.plot.axis)
            if (sum.plot.axis %in% c("p.value", "p.adj"))
                generateSummaryPlots(plot.data, file.name)
            else
                generateSummaryPlots(plot.data, file.name, 
                        Xlab = paste0("1 / ", sum.plot.axis))
        }
        file.name = paste0(summary.dir, sub(" - ", "-", 
                    contrast.names[i]), "-", gs.annot$label, "-methods")
        if (length(baseGSEAs) > 1 && !file.exists(paste0(file.name, 
                ".png")))
            generateMDSMethodsPlot(egsea.results[[i]], baseGSEAs, 
                    file.name)
    }
    
}

generateSumPlots.comparison <- function(egsea.results, egsea.comparison, gs.annot, 
gsa.dir,         
        sum.plot.cutoff=1, sum.plot.axis="p.value"){
    print("Comparison summary plots are being generated  ... ")
    
    summary.dir = paste0(gsa.dir, "/summary/")
    dir.create(file.path(summary.dir), showWarnings = FALSE)
    contrast.names = names(egsea.results)
    
    if (length(contrast.names) == 2){
        generateSummaryPlots.comparison(egsea.results, 
                egsea.comparison, gs.annot, 
                sum.plot.axis,file.prefix="", summary.dir)
    } 
    else if (length(contrast.names) > 2){
        print("There are more than 2 contrasts!")
        for (i in 1:(length(contrast.names)-1)){
            for (j in (i+1):length(contrast.names)){
                egsea.results.sub = egsea.results[c(i,j)]  
                generateSummaryPlots.comparison(egsea.results.sub, 
                        egsea.comparison, gs.annot, 
                        sum.plot.axis, file.prefix = 
                        paste0('-', i,j), summary.dir)
            }
        }
    }
}

generateSummaryPlots.comparison <- function(egsea.results, egsea.comparison, 
                        gs.annot, sum.plot.axis, file.prefix, summary.dir){
    file.name = paste0(summary.dir, gs.annot$label, file.prefix, "-summary")
    if (!file.exists(paste0(file.name, ".dir.png"))){
        contrast.names = names(egsea.results)
        plot.data = generatePlotData.comparison(egsea.results, 
                egsea.comparison, gs.annot, 
                sum.plot.axis)
        if (sum.plot.axis %in% c("p.value", "p.adj")){
            generateSummaryPlots(plot.data, file.name,
                    paste0("-log10(p-value) for ", 
                    contrast.names[1]),
                    paste0("-log10(p-value) for ", 
                    contrast.names[2]))
        }else{
            generateSummaryPlots(plot.data, file.name,
                    paste0("1 / ", sum.plot.axis, " for ", 
                        contrast.names[1]),
                    paste0("1 / ", sum.plot.axis, " for ", 
                        contrast.names[2]))
        }
        
    }
}

generatePlotData.comparison <- function(egsea.results.two, egsea.comparison, 
gs.annot,
        sum.plot.axis){ 
    
    gsets = as.character(rownames(egsea.comparison))    
#   print(length(egsea.results.all))
#   print(length(gsets))
#   print(head(fc))
    n = length(egsea.results.two)
    pvalues.all = matrix(0, length(gsets), n)
    gs.dirs.all = matrix(0, length(gsets), n)
    sig.combined = 0
    for (i in 1:n){
        egsea.results = egsea.results.two[[i]]      
        if ( sum.plot.axis %in% c("p.value", "p.adj")){
            pvalues = egsea.results[gsets, sum.plot.axis]       
    
            pvalues[pvalues == 0] = NA
            pvalues = -1 * log10(pvalues)
            pvalues[is.na(pvalues)] = max(pvalues, na.rm=TRUE) + 1  
    
        }else{
            pvalues = 1/egsea.results[gsets, sum.plot.axis] 
        }
        sig.combined = sig.combined + 
                egsea.results[gsets, "Significance"]
        pvalues.all[, i] = pvalues
        gs.dirs.all[, i] = egsea.results[gsets, "Direction"]
    }
    gs.dirs = rowMeans(gs.dirs.all, na.rm=TRUE)
    gs.sizes = as.numeric(sapply(
                    as.character(gs.annot$anno[gsets, "NumGenes"]), 
            function(x) strsplit(x, split="/", fixed=TRUE)[[1]][2]))
    sig.combined = sig.combined / n
    rank = seq(1, length(gs.dirs))
    
    plot.data = data.frame(id=gs.annot$anno[gsets, "ID"] , 
            x.data=pvalues.all[,1], y.data=pvalues.all[,2], 
            gsSize=gs.sizes, sig=sig.combined, dir=gs.dirs, rank = 
            rank)   
    plot.data = plot.data[order(plot.data[, "rank"]), ]
    return(plot.data)
}

generatePlotData <- function(egsea.results, gs.annot,
        sum.plot.cutoff,  sum.plot.axis){
    
    x.data = egsea.results[, sum.plot.axis]
    gsets = as.character(rownames(egsea.results))[x.data <= sum.plot.cutoff]
    x.data = x.data[x.data <= sum.plot.cutoff]
    if (sum.plot.axis %in% c("p.value", "p.adj")){
        x.data[x.data == 0] = NA
        x.data = -1 * log10(x.data)
        x.data[is.na(x.data)] = max(x.data, na.rm=TRUE) + 1     
    }
    else{
        x.data = 1 / x.data         
    }       
    
    rank = seq(1, length(gsets))

    gs.sizes = as.numeric(sapply(as.character(gs.annot$anno[gsets, 
"NumGenes"]), 
                    function(x) strsplit(x, split="/", 
fixed=TRUE)[[1]][2]))
    
    plot.data = data.frame(id=gs.annot$anno[gsets, "ID"] , 
            x.data=x.data, y.data=egsea.results[, "avg.logFC"], 
            gsSize=gs.sizes, sig=egsea.results[, "Significance"], 
            dir=egsea.results[, "Direction"], rank = rank)  
    return(plot.data)
}

# Generate MDS plot for the rankings of different methods
 
generateMDSMethodsPlot <- function(egsea.results, baseGSEAs, file.name, format=NULL){

    ranks = egsea.results[, baseGSEAs]
    if (is.null(format) || tolower(format) == "pdf"){
        pdf(paste0(file.name, ".pdf"), width=6, height=6)       
        limma::plotMDS(x=ranks, labels=baseGSEAs, col = "#4d4d4d", 
                xlab="Leading rank dim 1", 
                ylab="Leading rank dim 2", cex=1)
        dev.off()
    }
    if (is.null(format) || tolower(format) == "png"){
        png(paste0(file.name, ".png"), width=800, height=800)       
        limma::plotMDS(x=ranks, labels=baseGSEAs, col = "blue", 
                xlab="Leading rank dim 1", 
                ylab="Leading rank dim 2")
        dev.off()
    }
}

# Generate summary plots based on regulation direction and gene set ranking.
 
generateSummaryPlots <- function(plot.data, file.name, Xlab="-log10(p-value)",
        Ylab="Average Absolute logFC", format = NULL){     
    tryCatch({
        plot.data.sig = plot.data[plot.data[, "rank"] <= 10, ]
        sig.cols = rep("black", nrow(plot.data.sig))
    #   print(plot.data.sig)
    #       print(dim(plot.data))   
        # plot rank-based coloured bubbles
        p = qplot(x.data, y.data, data=plot.data, size=gsSize,asp=1, 
                colour=rank,
                xlab = Xlab, ylab = Ylab,
                xlim=c(0.8 * min(plot.data[, "x.data"], 
na.rm=TRUE), max(plot.data[, "x.data"], na.rm=TRUE)*1.05), 
                ylim=c(min(plot.data[, "y.data"], na.rm=TRUE), 
max(plot.data[, "y.data"], na.rm=TRUE) * 1.05))
        # customize bubbles colour 
        p = p + scale_colour_gradient(guide="colourbar", low="#56B1F7", 
high="#000000",
                limits=c(1,100), na.value="black", name="Rank")
        # customize bubble size
        p = p + scale_size("Cardinality", range=c(2,20))       
        if (is.null(format) || tolower(format) == "pdf"){
            pdf(paste0(file.name, ".rank.pdf"), width = 10, height = 7) 
        
            # label the bubbles of the top 10 gene sets
            print(p + geom_text(size=5, mapping=aes(x=x.data, y=y.data, 
                            label=id), 
                            data=plot.data.sig, 
                            colour=sig.cols, vjust=-1, hjust=1) )
            dev.off()       
        }
        if (is.null(format) || tolower(format) == "png"){
            png(paste0(file.name, ".rank.png"), width = 800, height = 700)
            print(p + geom_text(size=5, mapping=aes(x=x.data, y=y.data, 
    label=id), 
                            data=plot.data.sig, 
    colour=sig.cols, vjust=-1, hjust=1) )   
            dev.off()
        }
        # plot direction-based coloured bubbles
        top.10.ids = as.character(plot.data[plot.data[, "rank"] <= 10, 
"id"])
        sig.ids = setdiff(plot.data[rank(-plot.data[,"sig"], na.last = 
TRUE) <= 5, "id"], top.10.ids)
        sig.cols = c(rep("black", length(top.10.ids)), rep("blue", 
length(sig.ids)))
        plot.data.sig = plot.data[match(c(top.10.ids, sig.ids), 
plot.data[, "id"]), ]
        p = qplot(x.data, y.data, data=plot.data, size=sig,asp=1, 
                colour=dir,
                xlab = Xlab, ylab = Ylab,
                xlim=c(0.8 * min(plot.data[, "x.data"], 
na.rm=TRUE), max(plot.data[, "x.data"], na.rm=TRUE)*1.05),
                ylim=c(min(plot.data[, "y.data"], na.rm=TRUE), 
max(plot.data[, "y.data"], na.rm=TRUE) * 1.05))
        p = p + scale_colour_gradient(guide="colourbar", low="#56B1F7", 
high="#E35F5F",
                limits=c(-1,1), na.value="black", 
name="Regulation Direction") # low="#5FE377"
        p = p + scale_size("Significance", range=c(2,20))   
        if (is.null(format) || tolower(format) == "pdf"){
            pdf(paste0(file.name, ".dir.pdf"), width = 10, height = 7)  
            
            print(p + geom_text(size=5, mapping=aes(x=x.data, y=y.data, 
                            label=id), 
                            data=plot.data.sig, 
                            colour=sig.cols, vjust=-1, hjust=1) )
            dev.off()
        }
        if (is.null(format) || tolower(format) == "png"){
            png(paste0(file.name, ".dir.png"), width = 800, height = 700)
            print(p + geom_text(size=5, mapping=aes(x=x.data, y=y.data, 
                            label=id), 
                            data=plot.data.sig, 
                            colour=sig.cols, vjust=-1, hjust=1) )       
            dev.off()
        }
    }, 
    error = function(e){
        print(paste0("WARNING: summary plots were not generated for ", 
file.name))
    })
}
# Plot top KEGG pathways using the pathview package 
# 
# Plot top KEGG pathways using the pathview package  and overlay fold change 
# data on them.
 

plotPathways = function(gene.sets, fc, gs.annot, gsa.dir, kegg.dir, 
verbose=FALSE){
    #TODO: adjust the limit of the expression values 
    
    cat("Pathway maps are being generated for top-ranked \n pathways based 
			on logFC ... \n")
    current = getwd()   
    contrast.names = colnames(fc)
    if (is.null(kegg.dir))
        kegg.dir = paste0(gsa.dir, "/kegg-dir/")
    dir.create(file.path(kegg.dir))
    kegg.dir = normalizePath(kegg.dir)  

    for(i in 1:ncol(fc)){       
        print(paste0("    ", contrast.names[i]))
        pv.dir = paste0(gsa.dir, "/pv-top-gs-", gs.annot$label, "/", 
                sub(" - ", "-", contrast.names[i]))     
    
        dir.create(file.path(pv.dir), showWarnings = FALSE, recursive = 
            TRUE)
        setwd(pv.dir)
        for (gene.set in gene.sets[[i]]){
            id = as.character(gs.annot$anno[match(gene.set, 
                    gs.annot$anno[, "GeneSet"]), "ID"])       
            if (file.exists(paste0(pv.dir, "/", id, 
                ".pathview.png"))){                
                next
            }
            generatePathway(gene.set, gs.annot, fc[,i], kegg.dir, verbose)    
        }
    }
        
    setwd(current)
}

generatePathway <- function(gene.set, gs.annot, fc, kegg.dir="./", 
        verbose=FALSE, file.name = NULL){
    id = as.character(gs.annot$anno[match(gene.set, 
                            gs.annot$anno[, "GeneSet"]), "ID"]) 
    if (verbose)                
        pathview(gene.data=fc, pathway.id = id,  
                species = 
                        species.fullToShort[[tolower(gs.annot$species)]], 
                kegg.dir=kegg.dir, kegg.native=TRUE, 
                res=600,        
                limit=list(gene=c(-2,2), cpd=1), 
                bins=list(gene=20, cpd=10), 
                low = list(gene = "blue", cpd = 
                                "green"))
    else
        suppressMessages(
                pathview(gene.data=fc, pathway.id = id,  
                species = 
                        species.fullToShort[[tolower(gs.annot$species)]], 
                kegg.dir=kegg.dir, 
                kegg.native=TRUE, res=600,          
                limit=list(gene=c(-2,2), 
                        cpd=1), bins=list(gene=20, cpd=10), 
                low = list(gene = "blue", cpd = 
                                "green")))  
    if (!is.null(file.name)){
        pathway.file = paste0("./", id, ".pathview.png")
        if (file.exists(pathway.file)){
            file.rename(pathway.file, paste0(file.name, ".png"))
            file.remove(paste0(kegg.dir, id, ".png"))
            file.remove(paste0(kegg.dir, id, ".xml"))
        }else{
            stop("EGSEA could not generate the pathway map image.")
        }
    }
}

plotPathways.comparison <- function(gene.sets, fc, gs.annot, gsa.dir, kegg.dir, 
verbose=FALSE){
    cat("Pathway maps are being generated for top-ranked \n comparative 
			pathways based on logFC ... \n")
    current = getwd()
    if (is.null(kegg.dir))
        kegg.dir = paste0(gsa.dir, "/kegg-dir/")
    dir.create(file.path(kegg.dir))
    kegg.dir = normalizePath(kegg.dir)
    pv.dir = paste0(gsa.dir, "/pv-top-gs-", gs.annot$label, "/")    
    setwd(pv.dir)
    for (gene.set in gene.sets){
        id = gs.annot$anno[match(gene.set, gs.annot$anno[, "GeneSet"]), "ID"]
        if (file.exists(paste0(pv.dir, "/", id, ".pathview.multi.png"))){ 
            next
        }
        generateComparisonPathway(gene.set, gs.annot, fc, kegg.dir, verbose)
    }
    setwd(current)
}

generateComparisonPathway <- function(gene.set, gs.annot, fc, kegg.dir="./", 
        verbose=FALSE, file.name = NULL){
    id = as.character(gs.annot$anno[match(gene.set, 
                    gs.annot$anno[, "GeneSet"]), "ID"]) 
    if (!verbose)
        suppressMessages(pathview(gene.data=fc, pathway.id = id,  
                species = gs.annot$species, kegg.dir=kegg.dir, 
                kegg.native=TRUE, res=600, 
                limit=list(gene=c(-2,2), cpd=1), bins=list(gene=20, 
                        cpd=10), 
                key.align="y", key.pos="topright", 
                multi.state=TRUE, 
                low = list(gene = "blue", cpd = "green")))
    else
        pathview(gene.data=fc, pathway.id = id,  
                species = gs.annot$species, kegg.dir=kegg.dir, 
                kegg.native=TRUE, res=600, 
                limit=list(gene=c(-2,2), cpd=1), bins=list(gene=20, 
                        cpd=10), 
                key.align="y", key.pos="topright", multi.state=TRUE, 
                low = list(gene = "blue", cpd = "green"))
    
    if (!is.null(file.name)){
        pathway.file = paste0("./", id, ".pathview.multi.png")
        if (file.exists(pathway.file)){
            file.rename(pathway.file, paste0(file.name, ".png"))
            file.remove(paste0(kegg.dir, id, ".png"))
            file.remove(paste0(kegg.dir, id, ".xml"))
        }else{
            stop("EGSEA could not generate the pathway map image.")
        }
    }
}

plotHeatMapsLogFC.comparison <- function(gene.sets, fc, gs.annot,  symbolsMap, 
gsa.dir){   
    hm.dir = paste0(gsa.dir, "/hm-top-gs-", gs.annot$label, "/")        
    for (gene.set in gene.sets){
        id = gs.annot$anno[match(gene.set, gs.annot$anno[, "GeneSet"]), 
"ID"]
        file.name = paste0(hm.dir, "/", id, ".heatmap.multi.png")   
        #print(id)
        if (file.exists(file.name)){
            #print("  Heat map has been already generated.")
            next
        }           
        #print(head(fc))
        generateHeatMap(gene.set, gs.annot, fc, symbolsMap, sub(".png", 
"", file.name))             
    }
}

#Plot heat maps for each gene set in a gene sets annotation object 

plotHeatMapsLogFC = function(gene.sets, fc, gs.annot,  symbolsMap, gsa.dir){

    # create output directory for heatmaps
    print("Heat maps are being generated for top-ranked gene sets based on 
logFC ... ")
    if (!identical(rownames(fc), gs.annot$featureIDs)){     
        fc = fc[match(rownames(fc), gs.annot$featureIDs) , ]    
        if (!identical(rownames(fc), gs.annot$featureIDs)){
            print(rownames(fc)[1:10])       
            print(gs.annot$featureIDs[1:10])
            stop("The row names of the fold change matrix should 
match the featureIDs vector in the gs.annot list")          
        }
    }   
    if (!is.null(symbolsMap) && !identical(as.character(symbolsMap[,1]), 
gs.annot$featureIDs)){
        symbolsMap = symbolsMap[match(as.character(symbolsMap[,1]), 
gs.annot$featureIDs) , ]
        if (!identical(as.character(symbolsMap[,1]), 
gs.annot$featureIDs))
            stop("All featureIDs in the gs.annot list should map to 
a valid gene symbol")
    }
    
    contrast.names = colnames(fc)       
    
    for(i in 1:ncol(fc)){       
        hm.dir = paste0(gsa.dir, "/hm-top-gs-", gs.annot$label, "/", 
                sub(" - ", "-", contrast.names[i]))     
    
        dir.create(file.path(hm.dir), showWarnings = FALSE, recursive = 
TRUE)       
        
        for (gene.set in gene.sets[[i]]){
            id = gs.annot$anno[match(gene.set, gs.annot$anno[, 
"GeneSet"]), "ID"]
#           print(gene.set)
            file.name = paste0(hm.dir, "/", id, ".heatmap.png") 
        
            if (file.exists(file.name)){        
#               print("  Heat map has been already generated.")
                next
            }
            ##### generate heat map here    
            generateHeatMap(gene.set, gs.annot, fc[, i], 
symbolsMap, sub(".png", "", file.name))
        }
    }
        
}

generateHeatMap <- function(gene.set, gs.annot, fc, symbolsMap, file.name,
        format=NULL, print.csv=TRUE){    
    if (length(gs.annot$idx[[gene.set]]) < 2){
        warning(paste0("heatmap for ", gene.set, " is skipped. It has 
			only one gene."))
        return()
    }
    c = gs.annot$idx
    sel = which(names(c) == gene.set)
    sel.genes = c[[sel]]    
    if (!is.matrix(fc) || ncol(fc) == 1){
        hm = as.matrix(fc[sel.genes], nrow=length(fc[sel.genes]))
        hm = cbind(hm, hm)
        colnames(hm) = c(" ", " ")
    }else{
        hm = fc[sel.genes, ]
        colnames(hm) = gsub("X", "", colnames(fc))
    }
    if (is.null(symbolsMap))
        rownames(hm) = gs.annot$featureIDs[sel.genes] # genes in the 
# gene set    
    else
        rownames(hm) = symbolsMap[match(gs.annot$featureIDs[sel.genes], 
                                    symbolsMap[,1]), 2]
    
    
    # initialize and create the heat map plot
    colrange = colorpanel(74, "blue", "black", "red")
    br=c(seq(-2,-1,length=25),seq(-0.999,1,length=25), 
            seq(1.001,2,length=25))
    # adjust font size depending on the number of rows
    sel.genes.size = length(sel.genes)   
    if(sel.genes.size < 20){
        cr = 0.85
    }else if(sel.genes.size < 40){
        cr = 0.65
    }else if(sel.genes.size < 70){
        cr = 0.35
    }else if(sel.genes.size < 100){
        cr = 0.35
    }else 
        cr = 0.15
    # Export PDF file
    if (is.null(format) || tolower(format) == "pdf"){
        pdf(paste0(file.name, ".pdf"))
        par(cex.main = 0.55)
        heatmap.2(hm, col = colrange, breaks=br, 
                margins=c(10,10), cexRow=cr, 
                cexCol=0.85, trace = "none", 
                Colv = FALSE,  dendrogram="row", 
                main = paste0(gene.set, "(logFC)")
        ) 
        dev.off()   
    }
    # Export PNG file
    if (is.null(format) || tolower(format) == "png"){
        png(paste0(file.name, ".png")) # , width=800, height=800
        heatmap.2(hm, col = colrange, breaks=br, margins=c(10,10), 
                    cexRow=cr, cexCol=0.85, trace = "none", 
                    Colv = FALSE, Rowv=TRUE, dendrogram="none",
                    key=FALSE#, lmat=lmat, lwid = lwid, lhei = lhei 
                ) #,main = paste0(gene.set, " (logFC)"), colsep=c(3,6,9))   
        dev.off() 
    }
    # Export text file
    if (print.csv){
        if (!is.matrix(fc) || ncol(fc) == 1)
            write.csv(hm[,1], file=paste0(file.name, ".csv"), 
            row.names=TRUE, col.names=TRUE)
        else
            write.csv(hm, file=paste0(file.name, ".csv"), row.names=TRUE, 
            col.names=TRUE)
    }
}

#
#egsea.plotGenesetSizes <- function(gene.sets, gs.annot, fdr){
#   contrast.names = names(gene.sets)   
#   gs.idx = c()
#   labels = c()
#   cont.not.null = c()
#   
#   for (cont in contrast.names){
#       if (is.na(gene.sets[[cont]][1])){
#           cont.not.null = append(cont.not.null, FALSE)
#           next
#       }
#       if (length(gene.sets[[cont]]) > 15)
#           gsets = gene.sets[[cont]][1:15]
#       else
#           gsets = gene.sets[[cont]]
#       cont.not.null = append(cont.not.null, TRUE)
#       gs.idx = append(gs.idx, match(gsets, gs.annot$anno[,2]))    
    
#       labels = append(labels, rep(cont, length(gsets)))
#   }
#   if (length(labels) == 0){
#       print("The gene sets list is empty.")
#       return()
#   }
#       
#   sizes = as.numeric(gs.annot$anno[gs.idx, "NumGenes"]) 
#   gs.ids = gs.annot$anno[gs.idx, "ID"]
#   
#   no.classes = length(levels(factor(labels)))
#   if (no.classes > 4){
#       cols = rainbow(no.classes)
#   }else{
#       cols = c("tomato4", "seagreen4", "tomato1", 
# "seagreen2")[1:no.classes]
#   }
#   cols.labels = rep("black", length(gs.ids))
#   cols.labels[duplicated(gs.ids)] = "red" 
#   
#   rp = barplot(sizes, col=cols[factor(labels)],
#           xlim=c(0, max(sizes) * 1.2), ylab="Gene Sets", 
# xlab="Number of Genes",
#           names.arg=FALSE, horiz=TRUE,
#           main=paste0("Size of top-ranked gene sets in ", 
# gs.annot$label))
#   
#   if (length(gs.idx) < 30){
#       cex = 0.85
#   }else if (length(gs.idx) < 45){
#       cex = 0.65
#   }else 
#       cex = 0.5
#   
#   legend("topright", legend=paste(levels(factor(labels)), " (FDR=", 
# fdr[cont.not.null], ")"),
#           cex=0.5, fill=cols, text.col=cols)
#   text(0, rp, labels=gs.ids, cex=cex,
#           offset=1, adj=1.1, xpd=TRUE, col=cols.labels, pos=2)
#}
#
#egsea.plotGenesetRegDir <- function(gene.sets, test.results, gs.annot, fdr){
#   contrast.names = names(gene.sets)   
#   gs.idx = c()
#   labels = c()
#   dirs = c()  
#   cont.not.null = c()
#   pvalues = c()
#   for (cont in contrast.names){       
#       if (is.na(gene.sets[[cont]][1])){
#           cont.not.null = append(cont.not.null, FALSE)
#           next
#       }
#       if (length(gene.sets[[cont]]) > 15)
#           gsets = gene.sets[[cont]][1:15]
#       else
#           gsets = gene.sets[[cont]]
#       cont.not.null = append(cont.not.null, TRUE)
#       gs.idx = append(gs.idx, match(gsets, gs.annot$anno[,2]))    
    
#       labels = append(labels, rep(cont, length(gsets)))
#       dirs = append(dirs, test.results[[cont]][gsets, "Direction"])
#       pvalues = append(pvalues, test.results[[cont]][gsets, "PValue"])
#   }
#   if (length(labels) == 0){
#       print("The gene sets list is empty.")
#       return()
#   }
#   dirs = gsub("Up", +1, dirs)
#   dirs = gsub("Down", -1, dirs)
#   dirVal = as.numeric(dirs) # -1, +1
#   gs.ids = gs.annot$anno[gs.idx, "ID"]
#   
#   no.classes = length(levels(factor(labels)))
#   if (no.classes > 4){
#       cols = rainbow(no.classes)
#   }else{
#       cols = c("tomato4", "seagreen4", "tomato1", 
# "seagreen2")[1:no.classes]
#   }
#   
#   gs.ids.dup = gs.ids[duplicated(gs.ids)]
#   dirs.dup = dirVal[duplicated(gs.ids)]
#   labels.dup = labels[duplicated(gs.ids)]
#   pvalues.dup = pvalues[duplicated(gs.ids)]
#       
#   dirVal = dirVal[!(duplicated(gs.ids))]
#   labels = labels[!(duplicated(gs.ids))]
#   pvalues = pvalues[!duplicated(gs.ids)]
#   gs.ids = gs.ids[!(duplicated(gs.ids))]
#       
#   dirVal1 = rep(0, length(gs.ids))
#   dirVal1[match(gs.ids.dup, gs.ids)] = dirs.dup
#   labels1 = rep("", length(gs.ids))
#   labels1[match(gs.ids.dup, gs.ids)] = labels.dup
#   pvalues1 = rep(NA, length(gs.ids))
#   pvalues1[match(gs.ids.dup, gs.ids)] = pvalues.dup
#   
#   cols.labels = rep("black", length(gs.ids))
#   cols.labels[match(gs.ids.dup, gs.ids)] = "red"
#   
#   dirs = matrix(0, length(contrast.names), length(gs.ids))
#   pvals = matrix(NA, length(contrast.names), length(gs.ids))
#   colnames(dirs) = gs.ids
#   rownames(dirs) = contrast.names
#   colnames(pvals) = gs.ids
#   rownames(pvals) = contrast.names
#   for (i in 1:length(contrast.names)) {
#       dirs[i, which(labels == contrast.names[i])] = dirVal[labels == 
# contrast.names[i]]
#       dirs[i, which(labels1 == contrast.names[i])] = dirVal1 [labels1 
# == contrast.names[i]]
#       
#       pvals[i, which(labels == contrast.names[i])] = pvalues[labels 
# == contrast.names[i]]
#       pvals[i, which(labels1 == contrast.names[i])] = pvalues1 
# [labels1 == contrast.names[i]]
#   }
#   pvals = -1 * log2(pvals)
#   
#   rp = barplot(dirs[1, ], col=cols[1],
#           ylim=c(-1, 1.3), xlab="Gene Sets", ylab="Regulation 
#Direction",
#           names.arg=rep(NA, ncol(dirs)), beside=FALSE,
#           main=paste0("Regulation direction of top-ranked gene 
# sets in ", 
#                   gs.annot$label), add=FALSE)
#   rp = barplot(dirs[2, ], col=cols[2],
#           ylim=c(-1.3, 1.3), xlab="Gene Sets", ylab="Regulation 
# Direction",
#           names.arg=rep(NA, ncol(dirs)), beside=FALSE,
#           main=paste0("Regulation direction of top-ranked gene 
# sets in ", 
#                   gs.annot$label), add=TRUE)
#       
#   if (length(gs.idx) < 30){
#       cex = 0.85
#   }else if (length(gs.idx) < 45){
#       cex = 0.65
#   }else 
#       cex = 0.5
#   
#   legend("topright", legend=paste(levels(factor(labels)), " (FDR=", 
# fdr[cont.not.null], ")"), cex=0.8,
#           fill=cols, text.col=cols)
#   text(rp, par("usr")[3]+0.3, labels=gs.ids, cex=cex, srt=90,
#           offset=5,  pos=1, xpd=TRUE, col=cols.labels)
#   
#   
#   par(new=TRUE)
#   plot(1:ncol(pvals), pvals[1,],, 
# type="p",col="blue",xaxt="n",yaxt="n",xlab="",ylab="")
#   axis(4)
#   mtext("P-value",side=4,line=3)
#   par(new=TRUE)
#   plot(1:ncol(pvals), pvals[2,],, 
# type="p",col="red",xaxt="n",yaxt="n",xlab="",ylab="")
#   axis(4)
#   mtext("P-value",side=4,line=3)
#}


### Later the D3 Javascript library will be used to generated clickable plots 
