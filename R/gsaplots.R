#Ensemble of Gene Set Enrichment Analyses 
# 
# Author: Monther Alhamdoosh, E:m.hamdoosh@gmail.com
###############################################################################
# Plot GO graphs using topGO package

plotGOGraphs <- function(egsea.results, gs.annot, gsa.dir, sort.by){
    cat("GO graphs are being generated for top-ranked GO terms based on 
			p-values ... \n")        
    go.dir = paste0(gsa.dir, "/go-graphs/") 
    if (!dir.exists(go.dir))
        dir.create(file.path(go.dir), showWarnings = FALSE)
    contrast.names = names(egsea.results)
    file.name = paste0(go.dir, sub(" - ", "-", contrast.names), "-", 
        gs.annot@label, "-top-")
    if (file.exists(paste0(file.name[length(file.name)], "CC.png")))
        return()
    topGOdata = loadTopGOdata(gs.annot, go.dir)    
    noSig = 5 
    for (i in 1:length(contrast.names)){            
        print(contrast.names[i])
        generateGOGraphs(egsea.results[[i]], gs.annot, sort.by,
                file.name[i], topGOdata, noSig)
        
    }
    
}

plotGOGraphs.comparison <- function(egsea.results, gs.annot, gsa.dir, sort.by){
    cat("GO graphs are being generated for top-ranked COMPARISON\n 
		GO terms based on p-values ... \n")        
    go.dir = paste0(gsa.dir, "/go-graphs/") 
    if (!dir.exists(go.dir))
        dir.create(file.path(go.dir), showWarnings = FALSE)   
    file.name = paste0(go.dir, "comparison-", 
            gs.annot@label, "-top-")
    if (file.exists(paste0(file.name, "CC.png")))
        return()
    topGOdata = loadTopGOdata(gs.annot, go.dir)    
    noSig = 5
    generateGOGraphs(egsea.results, gs.annot, sort.by,
                file.name, topGOdata, noSig)        
   
}

loadTopGOdata <- function(gs.annot, go.dir=NULL){
    if (!is.null(go.dir) && file.exists(paste0(go.dir, "topGOdata.rda"))){
        load(paste0(go.dir, "topGOdata.rda"))
        return(topGOdata)
    }else{    
        if (tolower(gs.annot@species) %in% c("human", "homo sapiens")){
            mappingDB = "org.Hs.eg.db"
        }else if (tolower(gs.annot@species) %in% c("mouse", "mus musculus")){
            mappingDB = "org.Mm.eg.db"
        }else if (tolower(gs.annot@species) %in% c("rat", "rattus norvegicus")){
            mappingDB = "org.Rn.eg.db"
        }
        
        gl = rep(0, length(gs.annot@featureIDs))
        names(gl) = gs.annot@featureIDs
        capture.output(topGOdataBP <- new("topGOdata",ontology = "BP", allGenes = gl,
                        geneSel = topDiffGenes, nodeSize = 10, annot = 
                                annFUN.org, mapping=mappingDB))
        capture.output(topGOdataMF <- new("topGOdata",ontology = "MF", allGenes = gl,
                        geneSel = topDiffGenes, nodeSize = 10, annot = 
                                annFUN.org, mapping=mappingDB))
        capture.output(topGOdataCC <- new("topGOdata",ontology = "CC", allGenes = gl,
                        geneSel = topDiffGenes, nodeSize = 10, annot = 
                                annFUN.org, mapping=mappingDB)) 
        topGOdata = list(BP=topGOdataBP, MF=topGOdataMF, CC=topGOdataCC)
        if (!is.null(go.dir))
            save(topGOdata, file=paste0(go.dir, "topGOdata.rda"))
        return(topGOdata)
    }
}

#file.name : prefix
generateGOGraphs <- function(egsea.results, gs.annot, sort.by,
        file.name, topGOdata = NULL, noSig = 5, format = NULL){
    if (is.null(topGOdata)){
        topGOdata = loadTopGOdata(gs.annot)
    }
    go.subsets = list() 
    go.ids = as.character(gs.annot@anno[
                    match(rownames(egsea.results), 
                            gs.annot@anno[, "GeneSet"]), "GOID"])        
    go.subsets[["BP"]] = go.ids[go.ids %in% topGOdata[["BP"]]@graph@nodes]
    go.subsets[["MF"]] = go.ids[go.ids %in% topGOdata[["MF"]]@graph@nodes]
    go.subsets[["CC"]] = go.ids[go.ids %in% topGOdata[["CC"]]@graph@nodes]   
    
    scores = egsea.results[, sort.by]
    max.score = 1.001
    
    if (max(scores) > 1)
        scores = (scores - min(scores)) / (max(scores) - 
                    min(scores)) 
    scores = scores + 0.001
    names(scores) = as.character(
            gs.annot@anno[match(rownames(
                                    egsea.results), 
                            gs.annot@anno[, "GeneSet"]), 
                    "GOID"])       
    # Generate the BP graph
    tryCatch({
            scores.sub = rep(max.score, 
                    length(topGOdata[["BP"]]@graph@nodes))
            names(scores.sub) = topGOdata[["BP"]]@graph@nodes
            scores.sub[go.subsets[["BP"]]] = 
                    scores[go.subsets[["BP"]]]
            if (is.null(format) || tolower(format) == "pdf"){
                pdf(paste0(file.name, "BP.pdf"))     
                showSigOfNodes(topGOdata[["BP"]], scores.sub, 
                        firstSigNodes=noSig, 
                        useInfo='def', sigForAll=FALSE) # or 
                # use printGraph to write out plot to file
                dev.off()
            }
            if (is.null(format) || tolower(format) == "png"){
                png(paste0(file.name, "BP.png"), width=800, 
                        height=800)
                showSigOfNodes(topGOdata[["BP"]], scores.sub, 
                        firstSigNodes=noSig, 
                        useInfo='def', sigForAll=FALSE) # or 
                # use printGraph to write out plot to file
                dev.off()
            }
        }, error=function(err){
            print(paste("MY_ERROR:  ",err))
            dev.off()
            file.remove(paste0(file.name, "BP.pdf"))
            file.remove(paste0(file.name, "BP.png"))
        }
    )
#       stop("here")
    # write the MF graph
    tryCatch({
            scores.sub = rep(max.score, 
                    length(topGOdata[["MF"]]@graph@nodes))
            names(scores.sub) = topGOdata[["MF"]]@graph@nodes
            scores.sub[go.subsets[["MF"]]] = 
                    scores[go.subsets[["MF"]]]
            if (is.null(format) || tolower(format) == "pdf"){
                pdf(paste0(file.name, "MF.pdf"))
                showSigOfNodes(topGOdata[["MF"]], scores.sub, 
                        firstSigNodes=noSig, 
                        sigForAll=FALSE, useInfo='def') # or 
    # use printGraph to write out plot to file
                dev.off()
            }
            if (is.null(format) || tolower(format) == "png"){
                png(paste0(file.name, "MF.png"), width=800, 
                        height=800)
                showSigOfNodes(topGOdata[["MF"]], scores.sub, 
                        firstSigNodes=noSig, 
                        sigForAll=FALSE, useInfo='def') # or 
    # use printGraph to write out plot to file
                dev.off()
            }
        }, error=function(err){
            print(paste("MY_ERROR:  ",err))
            dev.off()
            file.remove(paste0(file.name, "MF.pdf"))
            file.remove(paste0(file.name, "MF.png"))
        }
    )       
    # write the CC graph
    tryCatch({
            scores.sub = rep(max.score, 
                    length(topGOdata[["CC"]]@graph@nodes))
            names(scores.sub) = topGOdata[["CC"]]@graph@nodes
            scores.sub[go.subsets[["CC"]]] = 
                    scores[go.subsets[["CC"]]]
            if (is.null(format) || tolower(format) == "pdf"){
                pdf(paste0(file.name, "CC.pdf"))
                showSigOfNodes(topGOdata[["CC"]], scores.sub, 
                        firstSigNodes=noSig, 
                        sigForAll=FALSE, useInfo='def') # or 
                # use printGraph to write out plot to file
                dev.off()
            }
            if (is.null(format) || tolower(format) == "png"){
                png(paste0(file.name, "CC.png"), width=800, 
                        height=800)
                showSigOfNodes(topGOdata[["CC"]], scores.sub, 
                        firstSigNodes=noSig, 
                        sigForAll=FALSE, useInfo='def') # or 
                # use printGraph to write out plot to file
                dev.off()
            }
        }, error=function(err){
            print(paste("MY_ERROR:  ",err))
            dev.off()
            file.remove(paste0(file.name, "CC.pdf"))
            file.remove(paste0(file.name, "CC.png"))
        }
    )
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
        contrast.names[i]), "-", gs.annot@label, "-summary")        

        if (!file.exists(paste0(file.name, ".dir.png"))){
            plot.data = generatePlotData(egsea.results[[i]], 
                     gs.annot, sum.plot.cutoff, sum.plot.axis)
            if (sum.plot.axis %in% c("p.value", "p.adj"))
                generateSummaryPlots(plot.data, file.name)
            else
                generateSummaryPlots(plot.data, file.name, 
                        Xlab = paste0("-", sum.plot.axis))
        }
        file.name = paste0(summary.dir, sub(" - ", "-", 
                    contrast.names[i]), "-", gs.annot@label, "-methods")
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
                sum.plot.axis, sum.plot.cutoff,
                file.prefix="", summary.dir)
    } 
    else if (length(contrast.names) > 2){
        print("There are more than 2 contrasts!")
        for (i in 1:(length(contrast.names)-1)){
            for (j in (i+1):length(contrast.names)){
                egsea.results.sub = egsea.results[c(i,j)]  
                generateSummaryPlots.comparison(egsea.results.sub, 
                        egsea.comparison, gs.annot, 
                        sum.plot.axis, sum.plot.cutoff, file.prefix = 
                        paste0('-', i,j), summary.dir)
            }
        }
    }
}

generateSummaryPlots.comparison <- function(egsea.results, egsea.comparison, 
                        gs.annot, sum.plot.axis, sum.plot.cutoff,
                        file.prefix, summary.dir){
    file.name = paste0(summary.dir, gs.annot@label, file.prefix, "-summary")
    if (!file.exists(paste0(file.name, ".dir.png"))){
        contrast.names = names(egsea.results)
        plot.data = generatePlotData.comparison(egsea.results, 
                egsea.comparison, gs.annot, 
                sum.plot.axis, sum.plot.cutoff)
        if (sum.plot.axis %in% c("p.value", "p.adj")){
            generateSummaryPlots(plot.data, file.name,
                    paste0("-log10(p-value) for ", 
                    contrast.names[1]),
                    paste0("-log10(p-value) for ", 
                    contrast.names[2]))
        }else{
            generateSummaryPlots(plot.data, file.name,
                    paste0("-", sum.plot.axis, " for ", 
                        contrast.names[1]),
                    paste0("-", sum.plot.axis, " for ", 
                        contrast.names[2]))
        }
        
    }
}

generatePlotData.comparison <- function(egsea.results.two, egsea.comparison, 
                 gs.annot,  sum.plot.axis, sum.plot.cutoff,
                 use.names = FALSE){ 
    tmp = egsea.results.two[[1]]
    gsets1 = as.character(rownames(tmp)) [tmp[, sum.plot.axis] <= sum.plot.cutoff]  
    tmp = egsea.results.two[[2]]
    gsets2 = as.character(rownames(tmp)) [tmp[, sum.plot.axis] <= sum.plot.cutoff]
    gsets = intersect(gsets1, gsets2)  
#    r1 = match(gsets, gsets1)
#    r2 = match(gsets, gsets2)
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
            pvalues[pvalues == 0] = 1*10^-22
            pvalues = -1 * log10(pvalues)
            pvalues[is.na(pvalues)] = max(pvalues, na.rm=TRUE) + 1 
        }else{
            pvalues =  - egsea.results[gsets, sum.plot.axis] 
        }
        sig.combined = sig.combined + 
                egsea.results[gsets, "Significance"]
        pvalues.all[, i] = pvalues
        gs.dirs.all[, i] = egsea.results[gsets, "Direction"]
    }
    gs.dirs = rowMeans(gs.dirs.all, na.rm=TRUE)
    gs.sizes = as.numeric(sapply(
                    as.character(gs.annot@anno[gsets, "NumGenes"]), 
            function(x) strsplit(x, split="/", fixed=TRUE)[[1]][2]))
    sig.combined = sig.combined / n
    rank = match(gsets, rownames(egsea.comparison))
    if (!use.names)
        labels = gs.annot@anno[gsets, "ID"] 
    else
        labels = gs.annot@anno[gsets, "GeneSet"]
    plot.data = data.frame(id=labels, 
            x.data=pvalues.all[,1], y.data=pvalues.all[,2], 
            gsSize=gs.sizes, sig=sig.combined, dir=gs.dirs, rank = 
            rank)   
    plot.data = plot.data[order(plot.data[, "rank"]), ]
    return(plot.data)
}

generatePlotData <- function(egsea.results, gs.annot,
        sum.plot.cutoff,  sum.plot.axis, use.names = FALSE){    
    x.data = egsea.results[, sum.plot.axis]
    gsets = as.character(rownames(egsea.results))[x.data <= sum.plot.cutoff]
    x.data = x.data[x.data <= sum.plot.cutoff]
    if (sum.plot.axis %in% c("p.value", "p.adj")){
        x.data[x.data == 0] = 1*10^-22
        x.data = -1 * log10(x.data)
        x.data[is.na(x.data)] = max(x.data, na.rm=TRUE) + 1     
    }
    else{
       x.data = - x.data         
    } 
    rank = seq(1, length(gsets))
    gs.sizes = as.numeric(sapply(as.character(gs.annot@anno[gsets, 
                "NumGenes"]), 
                function(x) strsplit(x, split="/", 
                    fixed=TRUE)[[1]][2]))
    if (!use.names)
        labels = gs.annot@anno[gsets, "ID"]
    else
        labels = gs.annot@anno[gsets, "GeneSet"]
    plot.data = data.frame(id=labels , 
            x.data=x.data, y.data=egsea.results[gsets, "avg.logFC"], 
            gsSize=gs.sizes, sig=egsea.results[gsets, "Significance"], 
            dir=egsea.results[gsets, "Direction"], rank = rank)  
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
        if (min(plot.data[, "x.data"], na.rm=TRUE) > 0){
            xlimits = c(0.8 * min(plot.data[, "x.data"], na.rm=TRUE), 
                max(plot.data[, "x.data"], na.rm=TRUE)*1.05)
        }else{
            xlimits = c(1.05 * min(plot.data[, "x.data"], na.rm=TRUE), 
                    max(plot.data[, "x.data"], na.rm=TRUE)*0.8)
        }
        if (max(plot.data[, "y.data"], na.rm=TRUE) > 0){
            ylimits = c(min(plot.data[, "y.data"], na.rm=TRUE), 
                    max(plot.data[, "y.data"], na.rm=TRUE) * 1.05)
        }else{
            ylimits = c(min(plot.data[, "y.data"], na.rm=TRUE), 
                    max(plot.data[, "y.data"], na.rm=TRUE) * 0.9)       
        }
    #   print(plot.data.sig)
    #       print(dim(plot.data))   
        # plot rank-based coloured bubbles
        p = qplot(x.data, y.data, data=plot.data, size=gsSize,asp=1, 
                colour=rank,
                xlab = Xlab, ylab = Ylab,
                xlim=xlimits, 
                ylim=ylimits)
        # customize bubbles colour 
        p = p + scale_colour_gradient(guide="colourbar", low="#56B1F7", 
high="#000000",
                limits=c(1,100), na.value="black", name="Rank")
        # customize bubble size
        p = p + scale_size("Cardinality", range=c(2,20))       
        if (is.null(format) || tolower(format) == "pdf"){
            pdf(paste0(file.name, ".rank.pdf"), width = 10, height = 7,
                    useDingbats = FALSE) 
        
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
                xlim=xlimits,
                ylim=ylimits)
        p = p + scale_colour_gradient(guide="colourbar", low="#56B1F7", 
high="#E35F5F",
                limits=c(-1,1), na.value="black", 
name="Regulation Direction") # low="#5FE377"
        p = p + scale_size("Significance", range=c(2,20))   
        if (is.null(format) || tolower(format) == "pdf"){
            pdf(paste0(file.name, ".dir.pdf"), width = 10, height = 7,
                    useDingbats = FALSE)  
            
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
    
    cat(paste0("Pathway maps are being generated for top-ranked \n pathways based ", 
"on logFC ... \n"))
    current = getwd()   
    contrast.names = colnames(fc)
    if (is.null(kegg.dir))
        kegg.dir = paste0(gsa.dir, "/kegg-dir/")
    dir.create(file.path(kegg.dir))
    kegg.dir = normalizePath(kegg.dir)  

    for(i in 1:ncol(fc)){       
        print(paste0("    ", contrast.names[i]))
        pv.dir = paste0(gsa.dir, "/pv-top-gs-", gs.annot@label, "/", 
                sub(" - ", "-", contrast.names[i]))     
    
        dir.create(file.path(pv.dir), showWarnings = FALSE, recursive = 
            TRUE)
        setwd(pv.dir)
        for (gene.set in gene.sets[[i]]){
            id = as.character(gs.annot@anno[match(gene.set, 
                    gs.annot@anno[, "GeneSet"]), "ID"])       
            if (file.exists(paste0(id,".pathview.png"))){                
                next
            }
            generatePathway(gene.set, gs.annot, fc[,i], kegg.dir, verbose)    
        }
    }
        
    setwd(current)
}

generatePathway <- function(gene.set, gs.annot, fc, kegg.dir="./", 
        verbose=FALSE, file.name = NULL){
    id = as.character(gs.annot@anno[match(gene.set, 
                            gs.annot@anno[, "GeneSet"]), "ID"]) 
    if (verbose)                
        pathview(gene.data=fc, pathway.id = id,  
                species = 
                        species.fullToShort[[tolower(gs.annot@species)]], 
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
                        species.fullToShort[[tolower(gs.annot@species)]], 
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
    cat(paste0("Pathway maps are being generated for top-ranked comparative\n", 
"pathways based on logFC ... \n"))
    current = getwd()
    if (is.null(kegg.dir))
        kegg.dir = paste0(gsa.dir, "/kegg-dir/")
    dir.create(file.path(kegg.dir))
    kegg.dir = normalizePath(kegg.dir)
    pv.dir = paste0(gsa.dir, "/pv-top-gs-", gs.annot@label, "/")    
    setwd(pv.dir)
    for (gene.set in gene.sets){
        id = gs.annot@anno[match(gene.set, gs.annot@anno[, "GeneSet"]), "ID"]
        if (file.exists(paste0(id, ".pathview.multi.png"))){ 
            next
        }
        generateComparisonPathway(gene.set, gs.annot, fc, kegg.dir, verbose)
    }
    setwd(current)
}

generateComparisonPathway <- function(gene.set, gs.annot, fc, kegg.dir="./", 
        verbose=FALSE, file.name = NULL){
    id = as.character(gs.annot@anno[match(gene.set, 
                    gs.annot@anno[, "GeneSet"]), "ID"]) 
    if (!verbose)
        suppressMessages(pathview(gene.data=fc, pathway.id = id,  
                species = gs.annot@species, kegg.dir=kegg.dir, 
                kegg.native=TRUE, res=600, 
                limit=list(gene=c(-2,2), cpd=1), bins=list(gene=20, 
                        cpd=10), 
                key.align="y", key.pos="topright", 
                multi.state=TRUE, 
                low = list(gene = "blue", cpd = "green")))
    else
        pathview(gene.data=fc, pathway.id = id,  
                species = gs.annot@species, kegg.dir=kegg.dir, 
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

plotHeatMapsLogFC.comparison <- function(gene.sets, fc, limma.tops, gs.annot,  symbolsMap, 
gsa.dir){   
    hm.dir = paste0(gsa.dir, "/hm-top-gs-", gs.annot@label, "/")        
    for (gene.set in gene.sets){
        id = gs.annot@anno[match(gene.set, gs.annot@anno[, "GeneSet"]), 
"ID"]
        file.name = paste0(hm.dir, "/", id, ".heatmap.multi.png")   
        #print(id)
        if (file.exists(file.name)){
            #print("  Heat map has been already generated.")
            next
        }           
        #print(head(fc))
        generateHeatMap(gene.set, gs.annot, fc, limma.tops, symbolsMap, sub(".png", 
"", file.name))             
    }
}

#Plot heat maps for each gene set in a gene sets annotation object 

plotHeatMapsLogFC = function(gene.sets, fc, limma.tops, gs.annot,  symbolsMap, gsa.dir){
    # create output directory for heatmaps
    print("Heat maps are being generated for top-ranked gene sets based on 
logFC ... ")
    if (!identical(rownames(fc), gs.annot@featureIDs)){     
        fc = fc[match(gs.annot@featureIDs, rownames(fc)) , ]    
        if (!identical(rownames(fc), gs.annot@featureIDs)){
            print(rownames(fc)[1:10])       
            print(gs.annot@featureIDs[1:10])
            stop("The row names of the fold change matrix should 
match the featureIDs vector in the gs.annot list")          
        }
    }   
    if (nrow(symbolsMap) > 0 && !identical(as.character(symbolsMap[,1]), 
                    gs.annot@featureIDs)){
        symbolsMap = symbolsMap[match(as.character(symbolsMap[,1]), 
                        gs.annot@featureIDs) , ]
        if (!identical(as.character(symbolsMap[,1]), 
                    gs.annot@featureIDs))
            stop("All featureIDs in the gs.annot list should map to 
					a valid gene symbol")
    }
    
    contrast.names = colnames(fc)       
    
    for(i in 1:ncol(fc)){       
        hm.dir = paste0(gsa.dir, "/hm-top-gs-", gs.annot@label, "/", 
            sub(" - ", "-", contrast.names[i]))
        dir.create(file.path(hm.dir), showWarnings = FALSE, 
            recursive = TRUE)       
        
        for (gene.set in gene.sets[[i]]){
            id = gs.annot@anno[match(gene.set, 
                gs.annot@anno[, "GeneSet"]), "ID"]
            file.name = paste0(hm.dir, "/", id, ".heatmap.png")
            if (file.exists(file.name)){        
                next
            }
            ##### generate heat map here    
            generateHeatMap(gene.set, gs.annot, fc[, i],
                limma.tops[contrast.names[i]],
                symbolsMap, sub(".png", "", file.name))
        }
    }
        
}

generateHeatMap <- function(gene.set, gs.annot, fc, limma.tops, symbolsMap, 
        file.name, format=NULL, print.csv=TRUE, 
        fc.colors= c("#67A9CF", "#F7F7F7", "#EF8A62")){
    stopifnot(length(fc.colors) == 3)
    if (length(gs.annot@idx[[gene.set]]) < 2){
        warning(paste0("heatmap for ", gene.set, " is skipped. It has 
			only one gene."))
        return()
    }
    idx = gs.annot@idx
    sel = which(names(idx) == gene.set)
    stopifnot(length(sel) == 1)
    sel.genes = idx[[sel]]    
    sig.genes = rep(FALSE, length(sel.genes))
    if (!is.matrix(fc) || ncol(fc) == 1){        
        tmpfc = fc[sel.genes]
        hm = as.matrix(tmpfc, nrow=length(tmpfc))
        hm = cbind(hm, hm)
        colnames(hm) = c(" ", " ")
        if (length(limma.tops) > 0){
            stopifnot(identical(names(fc), rownames(limma.tops[[1]])) ||
                    identical(rownames(fc), rownames(limma.tops[[1]])))
            csv.out = limma.tops[[1]][sel.genes, ]       
            sig.genes[csv.out[, "adj.P.Val"] <= 0.05] = TRUE
        }else
            csv.out = hm[,1]
        leglabels = c("FDR <= 0.05", 
                "FDR > 0.05")
        
    }else{
        hm = fc[sel.genes, ]        
#        colnames(hm) = gsub("X", "", colnames(fc))
        colnames(hm) = colnames(fc)
        if (length(limma.tops) > 0){
            csv.out = c()        
            for (contr in colnames(fc)){
                stopifnot(identical(rownames(fc), rownames(limma.tops[[contr]])))
                t = limma.tops[[contr]][sel.genes, ]
                sig.genes[t[, "adj.P.Val"] <= 0.05] = TRUE
                colnames(t) = gsub(" ", "", paste0(contr, ".", colnames(t)))
                if (length(csv.out) == 0)
                    csv.out = t
                else
                    csv.out = cbind(csv.out, t)                
            }
        }else{
            csv.out = hm
        }   
        leglabels = c("FDR <= 0.05 for at least one", 
                "FDR > 0.05 for all contrasts")
    }
    if (nrow(symbolsMap) == 0)
        rownames(hm) = gs.annot@featureIDs[sel.genes] # genes in the 
# gene set    
    else
        rownames(hm) = symbolsMap[match(gs.annot@featureIDs[sel.genes], 
                                    symbolsMap[,1]), 2]
                
    # initialize and create the heat map plot
    colrange = colorpanel(74, fc.colors[1], fc.colors[2],fc.colors[3])
    br=c(seq(-2,-0.5,length=25),seq(-0.499,0.5,length=25), 
            seq(0.51,2,length=25))
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
    # gene coloured based on significance 
    sig.color = rep("#6C6C6C", length(sel.genes))
    sig.color[sig.genes] = "#529F4E"
    # Export PDF file
    if (is.null(format) || tolower(format) == "pdf"){
        pdf(paste0(file.name, ".pdf"))
        if (nchar(gene.set) > 35)
            par(cex.main = 0.7)
        else
            par(cex.main = 0.8)
        heatmap.2(hm, col = colrange, breaks=br, 
                margins=c(10,10), cexRow=cr, 
                cexCol=0.85, trace = "none", 
                Colv = FALSE,  dendrogram="row", 
                density.info = "histogram",
                denscol = "black",
                #key.title = "",
                key.xlab = "logFC",
                key.ylab = "Count",
                colRow = sig.color,
                main = paste0(gene.set, "(logFC)")
        ) 
        legend("topright",      
                legend = leglabels,
                col = c("#529F4E", "#6C6C6C"), 
                title = "Significance of DE",
                lty= 1,             
                lwd = 5,           
                cex=.7
        )
        dev.off()   
    }
    # Export PNG file
    if (is.null(format) || tolower(format) == "png"){
        png(paste0(file.name, ".png")) # , width=800, height=800
        heatmap.2(hm, col = colrange, breaks=br, margins=c(10,10), 
                    cexRow=cr, cexCol=0.85, trace = "none", 
                    Colv = FALSE, Rowv=TRUE, dendrogram="none",
                    colRow = sig.color,
                    key=FALSE#, lmat=lmat, lwid = lwid, lhei = lhei 
                ) #,main = paste0(gene.set, " (logFC)"), colsep=c(3,6,9)) 
        dev.off() 
    }
    # Export text file
    if (print.csv){        
        write.csv(csv.out, file=paste0(file.name, ".csv"), 
            row.names=FALSE, col.names=TRUE)
    }
}


### Later the D3 Javascript library will be used to generated clickable plots 
