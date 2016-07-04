#Ensemble of Gene Set Enrichment Analyses 
# 
# Author: Monther Alhamdoosh, E:m.hamdoosh@gmail.com


egsea.main <- function(voom.results, contrast, gs.annots, baseGSEAs, 
                combineMethod, combineWeights,sort.by,  egsea.dir, 
                kegg.dir, logFC, symbolsMap, minSize, display.top, 
                logFC.cutoff, sum.plot.cutoff, sum.plot.axis, 
                vote.bin.width, print.base, verbose, num.threads, 
                report){
    baseGSEAs = sapply(baseGSEAs, tolower)
    baseGSEAs = unique(baseGSEAs)
    combineMethod = tolower(combineMethod)  
    sort.by = tolower(sort.by)
    if (sort.by == "significance")
        sort.by = "Significance"
    stopifnot(length(baseGSEAs) > 0 && length(setdiff(baseGSEAs, egsea.base())) == 0)
    stopifnot(combineMethod %in% egsea.combine())
    stopifnot(sort.by %in% c(egsea.sort()[1:8], baseGSEAs))
    if (length(baseGSEAs) <= 1){
        print("The ensemble mode was disabled.")
        if (sort.by %in% c("vote.rank", "avg.rank", "med.rank", 
               "min.pvalue", "min.rank")){
            sort.by = "p.value"
            warning("The sort.by argument was set to \"p.value\"")
        }
    }
    print("EGSEA analysis has started")
    # check arguments are valid
    if (!is.matrix(contrast)){
        stop("contrast argument must be a matrix object.")
    }
    if (is.null(colnames(contrast))){
        colnames(contrast) = paste0("contrast", rep(1, ncol(contrast)))
    }   
    # create output directory for 'egsea' results
    if (report){
        if (! dir.exists(file.path(egsea.dir))){
            dir.create(file.path(egsea.dir), showWarnings = FALSE)
            egsea.dir = normalizePath(egsea.dir)
        }
        ranked.gs.dir = paste0(egsea.dir, "/ranked-gene-sets-", combineMethod)
        dir.create(file.path(ranked.gs.dir), showWarnings = FALSE)  
        top.gs.dir = paste0(egsea.dir, "/top-gene-sets-", combineMethod)
        dir.create(file.path(top.gs.dir), showWarnings = FALSE)
    }
    else{
        top.gs.dir = NULL
        ranked.gs.dir = NULL
    }
        
    logFC.calculated = "No"
    if (is.null(logFC)){
        #row names should be Entrez Gene IDs in order to plot KEGG pathways
        logFC = getlogFCFromLMFit(voom.results, contrast)       
        logFC.calculated = "Yes"
    }else if (!is.matrix(logFC) || !identical(colnames(logFC), colnames(contrast))){
        stop("logFC should be a matrix object with column names equal 
			to the column names of the contrast matrix.")
    }
    
    if (class(gs.annots) == "GSCollectionIndex"){
        gs.annot = gs.annots
        gs.annots = list()
        gs.annots[[gs.annot@label]] = gs.annot
    }
    if (is.null(symbolsMap)){
        symbolsMap = data.frame()
    }else if (nrow(symbolsMap) > 0 && ncol(symbolsMap) >= 2){
        na.sym = is.na(symbolsMap[, 2])
        if (sum(na.sym) > 0){
            warning("Some \"NA\" Gene Symbols were replaced with Feature IDs")
            symbolsMap[na.sym, 2] = symbolsMap[na.sym, 1]
        }
    }
    gsas = EGSEAResults(contrasts = colnames(contrast), 
            sampleSize = getNumberofSamples(voom.results, contrast), 
            gs.annots = gs.annots, baseMethods=baseGSEAs,
            combineMethod = combineMethod, sort.by = sort.by,   
            symbolsMap = symbolsMap,
            logFC = logFC, report = report, report.dir = egsea.dir
            )
    
    skipped = c()
    for (gs.annot in gs.annots){
        gs.annot = selectGeneSets(gs.annot, min.size=minSize)  
        if (length(gs.annot@idx) == 0){
            print(paste0("No gene sets in ", gs.annot@label, " meets the 
minimum size criterion."))
            skipped = c(skipped, gs.annot@label)
            next
        }
        # run egsea and write out ranked gene sets for all contrasts
        print(paste0("EGSEA is running on the provided data and ",
                        gs.annot@label, " gene sets"))
        results.file = paste0(egsea.dir, "/", gs.annot@label,"-" , 
combineMethod , "-egsea-results.rda")
        if (file.exists(results.file)){     
            load(results.file)
            cat(paste0("The EGSEA results have been loaded from \n", 
results.file, "\n"))
            cat("If you want to re-run the EGSEA test,\n please remove this 
file or change the egsea.dir value.\n")
        }else{
            results <- runegsea(voom.results = voom.results, 
                    contrast=contrast, baseGSEAs=baseGSEAs,
                    combineMethod = combineMethod, combineWeights = NULL,
                    gs.annot = gs.annot, logFC = logFC, 
                    logFC.cutoff = logFC.cutoff, ranked.gs.dir = ranked.gs.dir, 
                    vote.bin.width=vote.bin.width, print.base = print.base,
                    report = report, num.workers = num.threads, verbose=verbose)            
            if (report){
                save(results, file=results.file)
            }
        }   
        egsea.results = results[["egsea.results"]]
        # order results based on the sort.by argument
        for (i in 1:length(egsea.results)){
            # sort based on the average ranking
            egsea.results[[i]] = egsea.results[[i]][
                                order(egsea.results[[i]][,sort.by],
                                decreasing=(sort.by == "Significance")), 
                                ]      
        }   
        
        # select top gene sets that pass an FDR cut-off threshold        
        gsets.top = egsea.selectTopGeneSets(egsea.results=egsea.results, 
                        fdr=display.top, gs.annot=gs.annot, 
                        top.gs.dir=top.gs.dir, report=report)   
        
        gsa = list("top.gene.sets"=gsets.top, "test.results"=egsea.results)
        if (print.base)
            gsa[["base.results"]] = results[["base.results"]]        
        
        # Generate heatmaps, pathways, GO graphs and summary plots
        if (report){
            if (length(grep("^kegg", gs.annot@label)) == 1){            
                plotPathways(gene.sets = gsets.top, fc=logFC, 
                        gs.annot=gs.annot, 
                        gsa.dir=egsea.dir, kegg.dir=kegg.dir, verbose = 
                        verbose)            
                file.name.pv = paste0(egsea.dir, "/pv-top-gs-", gs.annot@label, 
                        "/", 
                        sub(" - ", "-", colnames(contrast)),"-allPathways.html")
            }
            
            if ((gs.annot@label == "c5" || gs.annot@label == "gsdbgo") && 
                    "GOID" %in% colnames(gs.annot@anno)){
                plotGOGraphs(egsea.results=egsea.results,
                        gs.annot=gs.annot, gsa.dir=egsea.dir, sort.by=sort.by)
                file.name.go = paste0(egsea.dir, "/go-graphs/", 
                        sub(" - ", "-", colnames(contrast)),"-", gs.annot@label,
                        "-allGOgraphs.html")
            }
            
            plotHeatMapsLogFC(gene.sets = gsets.top, fc=logFC, 
                    gs.annot=gs.annot, 
                    symbolsMap=symbolsMap,
                    gsa.dir=egsea.dir)
            
            generateSumPlots(egsea.results = egsea.results, baseGSEAs = baseGSEAs, 
                    gs.annot = gs.annot, gsa.dir = egsea.dir,
                    sum.plot.cutoff = sum.plot.cutoff, sum.plot.axis = 
                    sum.plot.axis)
            
            # Select the annotations of the top gene sets and generate HTLM 
# reports 
            file.name = paste0(ranked.gs.dir, "/ranked-", gs.annot@label, 
                "-gene-sets-", 
                sub(" - ", "-", colnames(contrast)), '.txt')    
            file.name.hm = paste0(egsea.dir, "/hm-top-gs-", gs.annot@label, 
                "/", 
                sub(" - ", "-", colnames(contrast)),"-allHeatmaps.html")    
            file.name.sum = paste0(egsea.dir, "/summary/",  sub(" - ", "-", 
                colnames(contrast)),
                "-", gs.annot@label, "-summary.html")   
            
            # Create an HTML page for each contrast 
            for (i in 1:length(egsea.results)){         
                temp = egsea.results[[i]][1:ifelse(nrow(egsea.results[[i]]) > 
                    display.top, 
                    display.top, nrow(egsea.results[[i]])), ]
                writeEGSEAResultsToHTML(colnames(contrast)[i], temp, gs.annot, 
                    file.name[i])
                generateAllHeatmapsPage(colnames(contrast)[i], temp, gs.annot, 
                    file.name.hm[i])
                if (length(grep("^kegg", gs.annot@label)) == 1){
                    generateAllPathwaysPage(colnames(contrast)[i], temp, 
                        gs.annot, file.name.pv[i])
                }
                if ((gs.annot@label == "c5" || gs.annot@label == "gsdbgo") 
                        && "GOID" %in% colnames(gs.annot@anno) ){
                    generateAllGOgraphsPage(colnames(contrast)[i], gs.annot, 
                        file.name.go[i])
                }
                
                generateSummaryPage(colnames(contrast)[i], gs.annot, 
                    file.name.sum[i])
            }
        }
        
        
        # Comparison analysis reports generated here
        if (length(gsets.top) > 1){         
            egsea.comparison = createComparison(egsea.results, 
                    combineMethod = combineMethod, 
                    display.top=Inf, sort.by = sort.by)
            gsa[["comparison"]] = list()
            gsa$comparison[["test.results"]] = egsea.comparison
            egsea.comparison.all = egsea.comparison            
            egsea.comparison = egsea.comparison[1:ifelse(nrow(egsea.comparison) 
                > display.top, 
                display.top, nrow(egsea.comparison)), ]
            gsa$comparison[["top.gene.sets"]] = rownames(egsea.comparison)
            if (report){
                generateSumPlots.comparison(egsea.results = egsea.results, 
                    egsea.comparison = egsea.comparison.all, 
                    gs.annot = gs.annot, gsa.dir = egsea.dir, 
                    sum.plot.cutoff=sum.plot.cutoff, 
                    sum.plot.axis=sum.plot.axis)    
                plotHeatMapsLogFC.comparison(gene.sets = 
                    rownames(egsea.comparison), 
                    fc=logFC, gs.annot=gs.annot, symbolsMap=symbolsMap,
                    gsa.dir=egsea.dir)
                file.name = paste0(ranked.gs.dir, "/ranked-", gs.annot@label, 
                    "-gene-sets-compare.txt")   
                file.name.hm = paste0(egsea.dir, "/hm-top-gs-", gs.annot@label, 
                    "/allHeatmaps.html")                
                file.name.sum = paste0(egsea.dir, "/summary/",gs.annot@label, 
                    "-summary.html")                
                
                writeEGSEAResultsToHTML("Comparison Analysis", egsea.comparison
                        , gs.annot, file.name, comparison=TRUE)
                
                generateAllHeatmapsPage("Comparison Analysis", egsea.comparison
                        , gs.annot, file.name.hm, comparison=TRUE)
                
                if (length(grep("^kegg", gs.annot@label)) == 1){
                    plotPathways.comparison(gene.sets = 
                        rownames(egsea.comparison), fc=logFC, 
                        gs.annot=gs.annot, 
                        gsa.dir=egsea.dir, kegg.dir=kegg.dir, verbose = 
                        verbose)    
                    file.name.pv = paste0(egsea.dir, "/pv-top-gs-", 
                        gs.annot@label,"/allPathways.html") 
                    generateAllPathwaysPage("Comparison Analysis", 
                        egsea.comparison, 
                        gs.annot, file.name.pv, comparison=TRUE)
                }   
                if ((gs.annot@label == "c5" || gs.annot@label == "gsdbgo")  && 
                        "GOID" %in% colnames(gs.annot@anno)){
                    plotGOGraphs.comparison(egsea.results=egsea.comparison.all,
                            gs.annot=gs.annot, gsa.dir=egsea.dir, sort.by=sort.by)
                    file.name.go = paste0(egsea.dir, 
                        "/go-graphs/", gs.annot@label ,"-allGOgraphs.html")  
                    generateAllGOgraphsPage.comparison(colnames(contrast), 
                        gs.annot, file.name.go) 
                }
                
                generateSummaryPage.comparison(colnames(contrast), gs.annot, 
                    file.name.sum)      
            }
            
        }   
        gsas = addEGSEAResult(gsas, gs.annot@label, gsa)
        #gsas[[gs.annot@label]] = gsa
    }    
    if (report){
        gs.annots = gs.annots[! names(gs.annots) %in% skipped]
        generateEGSEAReport(voom.results, contrast, gs.annots, baseGSEAs, 
                            combineMethod, sort.by,  egsea.dir, kegg.dir, 
                            logFC.calculated, symbolsMap)    
        if (interactive()) try(browseURL(paste0("file://", 
                                normalizePath(egsea.dir),
                                    '/index.html')))
    }    
    print("EGSEA analysis has completed")
    return(gsas)
}


egsea.selectTopGeneSets <- function(egsea.results, fdr, gs.annot, 
top.gs.dir=NULL, report=TRUE){
    file.name = paste0(top.gs.dir, "/top-", gs.annot@label, "-gene-sets-", 
            sub(" - ", "-", names(egsea.results)), '.txt')
    contrast.names = names(egsea.results)
    gene.sets.fdr.detail = list()
    
    for(i in 1:length(egsea.results)){      
        #num.gene.sets.fdr = sum(egsea.results[[i]]$FDR < fdr[i])
        num.gene.sets.fdr = ifelse(length(gs.annot@idx) > fdr, 
fdr,length(gs.annot@idx))
        top.print = ifelse(length(gs.annot@idx) >= 10, 10,length(gs.annot@idx))
        if (top.print > num.gene.sets.fdr)
            top.print = num.gene.sets.fdr
        if (num.gene.sets.fdr > 0){
            egsea.results.top = egsea.results[[i]][1:num.gene.sets.fdr,]    
        
            gene.sets.fdr.detail[[i]] = rownames(egsea.results.top) 
            if (report){
                print(paste0("The top gene sets for contrast ", 
                                contrast.names[i], " are:"))
                if (length(grep("^kegg", gs.annot@label)) == 0){
                    top.table = cbind(gs.annot@anno[
                    match(rownames(egsea.results.top), gs.annot@anno[,2])
                                    ,-6],
                            egsea.results.top)
                    print(top.table[1:top.print, c("ID", "p.adj")])
                }else{
                    top.table = cbind(gs.annot@anno[
                                    match(rownames(egsea.results.top), 
    gs.annot@anno[,2])
                                    ,-2],
                            egsea.results.top)
                    print(top.table[1:top.print, c("Type", "p.adj")])
                }
            
                print(paste0("Writing out the top-ranked gene sets for each contrast .. 
                                        ", toupper(gs.annot@label), " gene sets"))
                write.table(top.table, file=file.name[i], sep="\t", quote=FALSE, 
row.names=FALSE)        
            }
            
        }else{
            write("No gene sets below current FDR found for this contrast, \n 
                            increase FDR or give up on this set." , 
file=file.name[i])
            print("No gene sets found")
            gene.sets.fdr.detail[[i]] = NA
        }
        
    }
    
    names(gene.sets.fdr.detail) = contrast.names
    return(gene.sets.fdr.detail)
}


runbaseGSEAParallelWorker <- function(args){
    #print(paste0("Running ", toupper(args$baseGSEA), " on all contrasts ... "))
    tryCatch({
                temp.result = runbaseGSEA(method=args$baseGSEA, 
                        args$voom.results, args$contrast, args$gs.annot,
                        args$ranked.gs.dir, output.base = args$print.base, 
                        num.threads = args$num.threads, verbose=args$verbose)
                if (args$verbose)
                    print(paste0("Running ", toupper(args$baseGSEA), " on all 
							contrasts ... COMPLETED "))
                else{
                    cat(args$baseGSEA)
                    cat("*")
                }
                return(temp.result)
            }, 
            error = function(e) {
                print(paste0("ERROR: ",toupper(args$baseGSEA), " encountered an 
							error ", e ))
            })
    return(NULL)
}

runegsea <- function(voom.results, contrast, baseGSEAs, combineMethod, 
                    combineWeights=NULL, gs.annot, logFC, logFC.cutoff, 
                    ranked.gs.dir, vote.bin.width, print.base=TRUE, 
                    report = TRUE, num.workers=8, verbose=FALSE){     
    # run egsea and write out ranked 'gene sets' for each 'contrast'    
    
    contrast.names = colnames(contrast)     
    # egsea.results.details stores Contrasts ==> Individual Results
    egsea.results.details = vector("list", ncol(contrast))  
    names(egsea.results.details) = contrast.names
    for (i in 1:length(contrast.names)){
        egsea.results.details[[i]] = vector("list", length(baseGSEAs))
        names(egsea.results.details[[i]]) = baseGSEAs
    }
    # run enrichment analysis using base methods
    args.all = list()
    threads.per.base = ceiling(num.workers / length(baseGSEAs))
    for (baseGSEA in baseGSEAs){
        args.all[[baseGSEA]] = list(baseGSEA=baseGSEA, 
                    voom.results=voom.results,
                    contrast=contrast, gs.annot=gs.annot, 
                    ranked.gs.dir=ranked.gs.dir,
                    print.base=(report && print.base), 
                    num.threads=threads.per.base,
                    verbose=verbose)
    }
    # temp.results stores Methods ==> Contrasts
    if (Sys.info()['sysname'] == "Windows" || num.workers <= 1 || 
length(baseGSEAs) == 1)
        # sequential processing
        temp.results = lapply(args.all, runbaseGSEAParallelWorker) 
    else
        # parallel processing
        temp.results = mclapply(args.all, runbaseGSEAParallelWorker, 
mc.cores=num.workers)   
    if (!verbose)
        cat("\n")
    # collect results   
#    print(baseGSEAs)
    for (baseGSEA in baseGSEAs){
        for (i in 1:ncol(contrast)){
            # order is important when combine
            if (is.null(temp.results[[baseGSEA]])){
                err = paste0("ERROR: One of the GSE methods failed on this 
dataset (",baseGSEA , ").\nRemove it and try again.\nSee error messages for 
more information.")
                stop(err)
            }
#            print(paste0(baseGSEA, colnames(contrast)[i]))
#            if (baseGSEA == "gage")
#                print(head(temp.results[[baseGSEA]][[i]]))
            egsea.results.details[[i]][[baseGSEA]] = 
temp.results[[baseGSEA]][[i]][names(gs.annot@idx),]         
        }
    }   
  
    # combine the results of base methods 
#    print(names(egsea.results.details))
    egsea.results = combineBaseGSEAs(results.multi=egsea.results.details, 
            combineMethod=combineMethod, combineWeights=combineWeights, 
            bin.width = vote.bin.width) 
    # calculate additional stats
#    print(names(egsea.results))
    for (i in 1:length(egsea.results)){ # i over contrasts
        gs.avg.fcs = numeric(0)     
        gs.dirs = numeric(0)
        gsets = as.character(rownames(egsea.results[[i]]))
        fc = logFC[, i]
        for (j in 1:length(gsets)){             
            sel.genes = gs.annot@idx[[gsets[j]]]
            gset.fc = fc[sel.genes]
            temp = abs(gset.fc)
            temp = temp[temp >= logFC.cutoff]
            gs.avg.fcs = c(gs.avg.fcs, mean(temp, na.rm=TRUE))          
            up = sum(gset.fc > logFC.cutoff, na.rm=TRUE)
            dn = sum(gset.fc < -logFC.cutoff, na.rm=TRUE)
            gs.dirs = c(gs.dirs, ifelse(up > dn, 1, -1))
        }   
#        print(head(egsea.results[[i]]))
        pvalues = egsea.results[[i]][, "p.adj"]
        pvalues[pvalues == 0] = 1*10^-22
        pvalues = -1 * log10(pvalues)
        pvalues[is.na(pvalues)] = max(pvalues, na.rm=TRUE) + 1  
        sig = pvalues * gs.avg.fcs
        if (max(sig, na.rm=TRUE) != min(sig, na.rm=TRUE))
            sig = (sig - min(sig, na.rm=TRUE)) / (max(sig, na.rm=TRUE) - 
                    min(sig, na.rm=TRUE)) * 100
        m = length(baseGSEAs)       
        if (m == 1){            
            egsea.results[[i]] = cbind(egsea.results[[i]], 
                        "avg.logFC"=gs.avg.fcs,
                    "Direction" = gs.dirs, "Significance" = sig)
        }
        else{
            n = ncol(egsea.results[[i]])
            # insert stat columns in the middle
            egsea.results[[i]] = cbind(egsea.results[[i]][, 1:(n-m)], 
                            "avg.logFC"=gs.avg.fcs,
                    "Direction" = gs.dirs, "Significance" = sig, 
                    egsea.results[[i]][, (n-m+1):n] )
        }
    }
    names(egsea.results) = colnames(contrast)
    results = list("egsea.results"=egsea.results)
    if (print.base){
        results[["base.results"]] = egsea.results.details
    }
    return(results)
}


combinePvalues <- function(data, combineMethod, combineWeights = NULL){
    if (ncol(data) > 1){
        if (combineMethod == "average"){
            #print(head(data))
            pvalues =sapply(1:nrow(data),  
                            function(i) {
                                x = data[i, ]                                
                                return(ifelse(length(x[!is.na(x)]) >= 4, 
                                        meanp(x[!is.na(x)])$p,
                                        mean(x)
                                      )) 
                             }
                           )    
        } else if (combineMethod == "fisher"){        
            data[data == 0] = 1*10^-22           
            pvalues = sapply(apply(data,  1, function(y) sumlog(y[!is.na(y)])), 
                    function(x) x$p)            
        } else if (combineMethod == "logitp"){        
            data[data == 0] = 1*10^-22
            data[data == 1] = 1 - 1*10^-5
            pvalues = sapply(apply(data,  1, function(y) logitp(y[!is.na(y)])), 
                    function(x) x$p)            
        }else if (combineMethod == "sump"){
            pvalues = sapply(apply(data,  1, function(y) sump(y[!is.na(y)])), 
                    function(x) x$p)            
        }else if (combineMethod == "sumz"){        
            data[data == 0] = 1*10^-22
            data[data == 1] = 1 - 1*10^-5
            pvalues = sapply(apply(data,  1, function(y) sumz(y[!is.na(y)])), #, combineWeights 
                    function(x) x$p)           
        }else if (combineMethod == "wilkinson"){        
            data[data == 0] = 1*10^-22
            data[data == 1] = 1 - 1*10^-5
            pvalues = sapply(apply(data,  1, function(y) wilkinsonp(y[!is.na(y)], r = 1)), 
                    function(x) x$p)            
        }    
    } else{
        pvalues = data[, 1]
    }
    #print(class(pvalues))
    #print(pvalues)
    adj.pvals = p.adjust(pvalues, method="BH")
    return(list(pvalues=pvalues, adj.pvals = adj.pvals))
}

combineBaseGSEAs <- function(results.multi, combineMethod, combineWeights=NULL, 
                    bin.width=5){
    if (length(results.multi[[1]]) == 1 && names(results.multi[[1]]) == "ora"){
        results.multi[[1]] = results.multi[[1]][[1]]
        results.multi[[1]] = results.multi[[1]][, 
                colnames(results.multi[[1]]) != "Rank"]
        return(results.multi)
    }
    # compress detailed results into matrices
    temp = extractPvaluesRanks(results.multi) 
    results.comp = temp$pvalues
    results.ranked = temp$ranks # gene sets ranked for each method
    results.combined = vector('list', length(results.ranked))
    names(results.combined) = names(results.comp)
#   results.comp = extractFDRs(results.multi)
    
    for (i in 1:length(results.comp)){ # i over contrasts
        temp = combinePvalues(results.comp[[i]], combineMethod, combineWeights)
        pvalues = temp$pvalues
        adj.pvals = temp$adj.pvals
        results.combined[[i]] = data.frame(cbind(
                "p.value"=pvalues, 
                "p.adj"=adj.pvals))
        if (ncol(results.ranked[[i]]) > 1){
            results.combined[[i]] = cbind(results.combined[[i]],  
                "vote.rank" = voteRank(results.ranked[[i]], 
                                        bin.width=bin.width),
                "avg.rank"=rowMeans(results.ranked[[i]], 
                            na.rm=TRUE), 
                "med.rank"=rowMedians(results.ranked[[i]], 
                            na.rm=TRUE),
                "min.pvalue"=sapply(1:nrow(results.comp[[i]]), 
                        function(x) min(results.comp[[i]][x, ], na.rm=TRUE)),
                "min.rank"=sapply(1:nrow(results.ranked[[i]]), 
                        function(x) min(results.ranked[[i]][x, ], na.rm=TRUE)),
                results.ranked[[i]]
            )
        }
        rownames(results.combined[[i]]) = rownames(results.comp[[i]])   
        
    }    
    
    #TODO: weighted average ranking
    #"wavg.rank"=rowMeans(sweep(results.ranked[[i]], MARGIN=2, combineWeights, 
    #`*`)), 
    
    return(results.combined)
}

voteRank <- function(results.ranked, bin.width=5){
    # convert to bins of bin.width 
    results.votes = numeric(nrow(results.ranked))
    if (bin.width != -1)
        results.ranked = (floor((results.ranked-1) / bin.width) + 1) * bin.width
    
    # apply majority voting: find frequencies. if tie, ignore 
    for (i in 1:nrow(results.ranked)){
        votes = results.ranked[i, ]
        cnts = list()
        for (j in 1:length(votes)){
            if (!is.null(cnts[[paste0(votes[j])]]))
                cnts[[paste0(votes[j])]] = cnts[[paste0(votes[j])]] + 1
            else
                cnts[[paste0(votes[j])]] = 1
        }       
        #if (max(as.numeric(cnts)) > floor(ncol(results.ranked) / 2)){
        indx = which.max(as.numeric(cnts))
        results.votes[i] = names(cnts)[indx]
#       }else
#           results.votes[i] = NA
    }
    # return the new ranking
    return(as.numeric(results.votes))
}



extractAdjPvaluesRanks <- function(results.multi){
    results.comp = vector("list", length(results.multi))
    names(results.comp) = names(results.multi)
    results.ranked = vector("list", length(results.comp))
    names(results.ranked) = names(results.comp)
    for (i in 1:length(results.multi)){ # i over contrasts
        results.comp[[i]] = matrix(0, nrow(results.multi[[i]][[1]]), 
                length(results.multi[[i]]))
        rownames(results.comp[[i]]) = rownames(results.multi[[i]][[1]])
        colnames(results.comp[[i]]) = names(results.multi[[i]]) 
        results.ranked[[i]] = matrix(0, nrow(results.comp[[i]]), 
                ncol(results.comp[[i]]))    
        rownames(results.ranked[[i]]) = rownames(results.comp[[i]])
        colnames(results.ranked[[i]]) = colnames(results.comp[[i]])
        j = 1
        #print(names(results.comp)[i])
        for (method in names(results.multi[[i]])){ # j is over base methods
            #print(method)
            #print(summary(results.multi[[i]][[j]]))
            if (method %in% c("camera", "roast", "fry")){
                results.comp[[i]][,j] = as.numeric(results.multi[[i]][[j]][, 
                                "FDR"])  
            }else if (method == "gage"){
                results.comp[[i]][,j] = as.numeric(results.multi[[i]][[j]][, 
                                "q.val"])
            }else if (method == "padog"){
                results.comp[[i]][,j] = as.numeric(results.multi[[i]][[j]][, 
                                "Ppadog"])
            }else if (method %in% c("plage","zscore", "gsva", "ssgsea")){
                results.comp[[i]][,j] = as.numeric(results.multi[[i]][[j]][, 
                                "adj.P.Val"]) 
            }else if (method == "globaltest"){
                results.comp[[i]][,j] = as.numeric(results.multi[[i]][[j]][, 
                                "p-value"])
            }else if (method == "ora"){
                results.comp[[i]][,j] = as.numeric(results.multi[[i]][[j]][, 
                                "p.adj"])
            }else if (method == "safe"){
                results.comp[[i]][,j] = as.numeric(results.multi[[i]][[j]][, 
                                "Adj.p.value"]) 
            }else if (method == "spia"){
                results.comp[[i]][,j] = as.numeric(results.multi[[i]][[j]][, 
                                "pG"])
            }
            results.ranked[[i]][,j] = as.numeric(results.multi[[i]][[j]][, 
                            "Rank"])
            if (max(results.comp[[i]][,j], na.rm = TRUE) < 0.000001)
                print(paste0("WARNING: ", method, " produces very low p-values 
                                        on ", names(results.comp)[i]))
            j = j + 1
        }
    }
    return(list(pvalues=results.comp, ranks=results.ranked))
}

extractPvaluesRanks <- function(results.multi){
    results.comp = vector("list", length(results.multi))
    names(results.comp) = names(results.multi)
    results.ranked = vector("list", length(results.comp))
    names(results.ranked) = names(results.comp)
    for (i in 1:length(results.multi)){ # i over contrasts
        results.comp[[i]] = matrix(0, nrow(results.multi[[i]][[1]]), 
length(results.multi[[i]]))
        rownames(results.comp[[i]]) = rownames(results.multi[[i]][[1]])
        colnames(results.comp[[i]]) = names(results.multi[[i]]) 
        results.ranked[[i]] = matrix(0, nrow(results.comp[[i]]), 
ncol(results.comp[[i]]))    
        rownames(results.ranked[[i]]) = rownames(results.comp[[i]])
        colnames(results.ranked[[i]]) = colnames(results.comp[[i]])
        j = 1
        #print(names(results.comp)[i])
        for (method in names(results.multi[[i]])){ # j is over base methods
            #print(method)
            #print(summary(results.multi[[i]][[j]]))
            if (method %in% c("camera", "roast", "fry")){
                results.comp[[i]][,j] = as.numeric(results.multi[[i]][[j]][, 
"PValue"])  
            }else if (method == "gage"){
                results.comp[[i]][,j] = as.numeric(results.multi[[i]][[j]][, 
"p.val"])
            }else if (method == "padog"){
                results.comp[[i]][,j] = as.numeric(results.multi[[i]][[j]][, 
"Ppadog"])
            }else if (method %in% c("plage","zscore", "gsva", "ssgsea")){
                results.comp[[i]][,j] = as.numeric(results.multi[[i]][[j]][, 
"P.Value"]) 
            }else if (method == "globaltest"){
                results.comp[[i]][,j] = as.numeric(results.multi[[i]][[j]][, 
"p-value"])
            }else if (method == "ora"){
                results.comp[[i]][,j] = as.numeric(results.multi[[i]][[j]][, 
"p.value"])
            }else if (method == "safe"){
                results.comp[[i]][,j] = as.numeric(results.multi[[i]][[j]][, 
"P.value"]) 
            }else if (method == "spia"){
                results.comp[[i]][,j] = as.numeric(results.multi[[i]][[j]][, 
"pG"])
            }
            results.ranked[[i]][,j] = as.numeric(results.multi[[i]][[j]][, 
"Rank"])
            if (max(results.comp[[i]][,j], na.rm = TRUE) < 0.000001)
                print(paste0("WARNING: ", method, " produces very low p-values 
on ", names(results.comp)[i]))
            j = j + 1
        }
    }
    return(list(pvalues=results.comp, ranks=results.ranked))
}

runbaseGSEA <- function(method, voom.results, contrast, gs.annot, 
                        ranked.gs.dir, output.base=TRUE,
                        num.threads = 4, verbose=TRUE){
    if (method == "camera"){        
        return(runcamera(voom.results = voom.results, contrast = contrast, 
gs.annot = gs.annot, 
                        ranked.gs.dir= ranked.gs.dir, output = output.base, 
num.workers = num.threads, verbose = verbose))     
    }else if (method == "roast"){       
        return(runroast(voom.results = voom.results, contrast = contrast, 
gs.annot = gs.annot, 
                        ranked.gs.dir= ranked.gs.dir, output = output.base, 
num.workers = num.threads, verbose = verbose))     
    }else if (method == "fry"){
        return(runfry(voom.results = voom.results, contrast = contrast, 
                        gs.annot = gs.annot, 
                        ranked.gs.dir= ranked.gs.dir, output = output.base, 
                        num.workers = num.threads, verbose = verbose))
    }else if (method == "gage"){        
        return(rungage(voom.results = voom.results, contrast = contrast, 
gs.annot = gs.annot, 
                        ranked.gs.dir= ranked.gs.dir, output = output.base, 
num.workers = num.threads, verbose = verbose))     
    }else if (method == "padog"){       
        return(runpadog(voom.results = voom.results, contrast = contrast, 
gs.annot = gs.annot, 
                        ranked.gs.dir= ranked.gs.dir, output = output.base, 
num.workers = num.threads, verbose = verbose))
    }else if (method == "plage"){       
        return(rungsva(method="plage", voom.results = voom.results, contrast = 
contrast, gs.annot = gs.annot, 
                        ranked.gs.dir= ranked.gs.dir, output = output.base, 
num.workers = num.threads, verbose = verbose))
    }else if (method == "zscore"){  
        return(rungsva(method="zscore", voom.results = voom.results, contrast = 
contrast, gs.annot = gs.annot, 
                        ranked.gs.dir= ranked.gs.dir, output = output.base, 
num.workers = num.threads, verbose = verbose))
    }else if (method == "gsva"){        
        return(rungsva(method="gsva", voom.results = voom.results, contrast = 
contrast, gs.annot = gs.annot, 
                        ranked.gs.dir= ranked.gs.dir, output = output.base, 
num.workers = num.threads, verbose = verbose))
    }else if (method == "ssgsea"){  
        return(rungsva(method="ssgsea", voom.results = voom.results, contrast = 
contrast, gs.annot = gs.annot, 
                        ranked.gs.dir= ranked.gs.dir, output = output.base, 
num.workers = num.threads, verbose = verbose))
    }else if (method == "globaltest"){      
        return(runglobaltest(voom.results = voom.results, contrast = contrast, 
gs.annot = gs.annot, 
                        ranked.gs.dir= ranked.gs.dir, output = output.base, 
num.workers = num.threads, verbose = verbose))
    }else if (method == "safe"){        
        return(runsafe(voom.results = voom.results, contrast = contrast, 
gs.annot = gs.annot, 
                        ranked.gs.dir= ranked.gs.dir, output = output.base, 
num.workers = num.threads, verbose = verbose))
    }else if (method == "ora"){
        return (runora(voom.results = voom.results, contrast = contrast, 
gs.annot = gs.annot, 
                        ranked.gs.dir= ranked.gs.dir, output = output.base, 
num.workers = num.threads, verbose = verbose))         
    }else{
        stop("Method not recognized. Type egsea.base() to see supported base 
methods.")
    }
}




createComparison <- function(egsea.results, combineMethod="fisher", display.top=100, 
        sort.by="p.value"){
    egsea.comparison = numeric(0)
    col.names = colnames(egsea.results[[1]])
    gset.names =  rownames(egsea.results[[1]])
    for (i in 1:length(col.names)){ # iterate over egsea.results columns
        if (col.names[i] == "p.adj")
            next
        temp= numeric(0)
        for (j in 1:length(egsea.results)){ # iterate over contrasts
            temp = cbind(temp, egsea.results[[j]][gset.names, i])
        }
        
        if (col.names[i] == "p.value"){         
            temp = combinePvalues(temp, combineMethod)
            pvalues = temp$pvalues
            adj.pvals = temp$adj.pvals
            egsea.comparison = cbind(egsea.comparison, 
                    pvalues, adj.pvals)
        }
        else if (col.names[i] == "med.rank"){
            egsea.comparison = cbind(egsea.comparison, rowMedians(temp, 
                    na.rm=TRUE))
        }
        else if (length(grep("min", col.names[i])) > 0){
            minVals = sapply(1:nrow(temp), function(x) min(temp[x, ], 
                    na.rm=TRUE))
            egsea.comparison = cbind(egsea.comparison, minVals)
        }
        else if (col.names[i] %in% c("avg.rank", "Direction", "Significance", 
                "avg.logFC"))
            egsea.comparison = cbind(egsea.comparison, rowMeans(temp, 
                na.rm=TRUE))
        else if (col.names[i] == "vote.rank"){
            if (length(egsea.results) > 2)
                egsea.comparison = cbind(egsea.comparison, voteRank(temp, 
                    bin.width = -1))
            else
                egsea.comparison = cbind(egsea.comparison, rowMeans(temp, 
                            na.rm=TRUE))
        }
    }
    rownames(egsea.comparison) = gset.names
#   print(colnames(egsea.comparison))
#   print(col.names)
    colnames(egsea.comparison) = col.names[1:length(colnames(egsea.comparison))]
    egsea.comparison = egsea.comparison[order(egsea.comparison[,sort.by]), ]
    display.top = ifelse(nrow(egsea.comparison) > display.top, display.top, 
nrow(egsea.comparison))
    return(egsea.comparison[1:display.top, ])
    
}

getlogFCFromLMFit <- function(voom.results, contrast){
    # to be changed for gene symbols support
    print("Log fold changes are estimated using limma package ... ")
    # fit linear model for each gene using limma package functions
    vfit = lmFit(voom.results, design=voom.results$design) # Fit linear model 
# for each gene given a series of arrays
    vfit = contrasts.fit(vfit, contrast) # make all pair-wise comparisons 
# between the groups
    ebayes.results = eBayes(vfit) # compute moderated t-statistics, moderated 
#F-statistic, 
    #and log-odds of differential expression by empirical Bayes moderation of 
# the 
    #standard errors towards a common value         
    
    logFC = matrix(0, nrow(ebayes.results), ncol(contrast))
    for (i in 1:ncol(contrast)){
        top.table = topTable(ebayes.results, coef=i, number=Inf, sort.by="none")
        limma.fc = top.table$logFC      
        names(limma.fc) = rownames(ebayes.results)  
        logFC[, i] = limma.fc
    }
    
    rownames(logFC) = rownames(ebayes.results)
    colnames(logFC) = colnames(contrast)
    return(logFC)
}
