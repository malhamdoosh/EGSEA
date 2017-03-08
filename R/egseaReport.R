# TODO: Add comment
# 
# Author: monther
###############################################################################


generate.EGSEA.Report <- function(egseaResults, limma.tops = list(), 
        display.top = 20, sort.by = NULL, 
        egsea.dir = NULL, kegg.dir = NULL, 
        sum.plot.axis = NULL, sum.plot.cutoff = NULL, 
        num.threads = 4,
        print.base = FALSE,
        verbose = FALSE){
    # create output directory for 'egsea' results
    print("EGSEA HTML report is being generated ...")
    start.time <- proc.time()
    timestamp()
    if (is.null(egsea.dir))
        egsea.dir = egseaResults@report.dir
    if (is.null(sort.by))
        sort.by = egseaResults@sort.by
    if (is.null(sum.plot.axis))
        sum.plot.axis = egseaResults@sum.plot.axis
    if (is.null(sum.plot.cutoff))
        sum.plot.cutoff = egseaResults@sum.plot.cutoff  
    contr.names = egseaResults@contr.names
    if (length(limma.tops) == 0 && 
            length(egseaResults@limmaResults) > 0){
        limma.tops = get.toptables(egseaResults@limmaResults,
                egseaResults@contrast)       
    }
    combineMethod = egseaResults@combineMethod
    baseGSEAs = egseaResults@baseMethods
    logFC = egseaResults@logFC
    logFC.calculated = egseaResults@logFC.calculated
    gs.annots = egseaResults@gs.annots
    symbolsMap = egseaResults@symbolsMap
    
    if (is.null(egsea.dir)){
        timestamp = as.integer(as.numeric(Sys.time()))
        egsea.dir = paste0("./egsea_results_", timestamp)
    }
    if (! dir.exists(file.path(egsea.dir))){
        dir.create(file.path(egsea.dir), showWarnings = FALSE)        
    }else{
        cat("WARNING: EGSEA report has been found. Change 'egsea.dir' to generate a new report.\n")
    }
    egsea.dir = normalizePath(egsea.dir)
    ranked.gs.dir = paste0(egsea.dir, "/ranked-gene-sets-base")
    dir.create(file.path(ranked.gs.dir), showWarnings = FALSE)
    summary.dir = paste0(egsea.dir, "/summary/")
    dir.create(file.path(summary.dir), showWarnings = FALSE)
    gs.annots = gs.annots[names(gs.annots) %in% names(egseaResults@results)]
    for (gs.annot in gs.annots){
        print(paste0("Report pages and figures are being generated for the ",
                        gs.annot@label, " collection ...")) 
        egsea.results = egseaResults@results[[gs.annot@label]]$test.results
        base.results = egseaResults@results[[gs.annot@label]]$base.results
        # order results based on the sort.by argument
        for (i in 1:length(egsea.results)){
            # sort based on the average ranking
            egsea.results[[i]] = egsea.results[[i]][
                    order(egsea.results[[i]][,sort.by],
                            decreasing=(sort.by == "significance")), 
            ]      
        }  
        # select top gene sets   
        gsets.top = egsea.selectTopGeneSets(egsea.results = egsea.results, 
                display.top = display.top, gs.annot = gs.annot, 
                report=TRUE) 
        # write out the base statistics if available 
        if (!is.null(base.results) && print.base){            
            write.base.stats(base.results, baseGSEAs, contr.names,
                    gs.annot, ranked.gs.dir, num.threads)
        }
        # generate summary heatmap for each collection 
        file.name.sum = paste0(summary.dir, gs.annot@label, 
                "-summary-heatmap-", sort.by) 
        if (!file.exists(paste0(file.name.sum, ".png")))
            plotSummaryHeatmap(egseaResults, gs.annot$label, 
                number = display.top,
                sort.by = sort.by, hm.vals = "avg.logfc.dir",  
                file.name = file.name.sum, 
                format = NULL, verbose = verbose)
        # Generate heatmaps, pathways, GO graphs and summary plots
        if (length(grep("^kegg", gs.annot@label)) == 1){            
            plotPathways(gene.sets = gsets.top, fc=logFC, 
                    gs.annot=gs.annot, 
                    gsa.dir=egsea.dir, kegg.dir=kegg.dir, verbose = 
                            verbose)            
            file.name.pv = paste0(egsea.dir, "/pv-top-gs-", gs.annot@label, 
                    "/", 
                    contr.names,"-allPathways.html")
        }
        
        if ((gs.annot@label == "c5" || gs.annot@label == "gsdbgo") && 
                "GOID" %in% colnames(gs.annot@anno)){
            plotGOGraphs(egsea.results=egsea.results,
                    gs.annot=gs.annot, gsa.dir=egsea.dir, sort.by=sort.by,
                    verbose)
            file.name.go = paste0(egsea.dir, "/go-graphs/", 
                    contr.names,"-", gs.annot@label,
                    "-allGOgraphs.html")
        }
        
        plotHeatMapsLogFC(gene.sets = gsets.top, fc=logFC, 
                limma.tops = limma.tops,
                gs.annot=gs.annot, 
                symbolsMap=symbolsMap,
                gsa.dir=egsea.dir)
        
        generateSumPlots(egsea.results = egsea.results, baseGSEAs = baseGSEAs, 
                gs.annot = gs.annot, gsa.dir = egsea.dir,
                sum.plot.cutoff = sum.plot.cutoff, 
                sum.plot.axis = sum.plot.axis)
        
        # Select the annotations of the top gene sets and generate HTLM 
# reports 
        file.name = paste0(ranked.gs.dir, "/ranked-", gs.annot@label, 
                "-gene-sets-", 
                contr.names, '.txt')    
        file.name.hm = paste0(egsea.dir, "/hm-top-gs-", gs.annot@label, 
                "/", 
                contr.names,"-allHeatmaps.html")    
        file.name.sum = paste0(summary.dir,  sub(" - ", "-", 
                        contr.names),
                "-", gs.annot@label, "-summary.html")   
        
        # Create an HTML page for each contrast 
        for (i in 1:length(egsea.results)){  
            file.name.bar = paste0(summary.dir,                      
                    contr.names [i], "-", gs.annot@label, "-bar-plot-",
                    sort.by) 
            if (!file.exists(paste0(file.name.bar, ".png")))
                plotBars(egseaResults, gs.annot@label, contrast = i, 
                    number = display.top, sort.by = sort.by,
                    bar.vals = "p.adj",
                    file.name = file.name.bar, 
                    format = NULL, verbose=verbose)
            temp = egsea.results[[i]][1:ifelse(nrow(egsea.results[[i]]) > 
                                    display.top, 
                            display.top, nrow(egsea.results[[i]])), ]
            writeEGSEAResultsToHTML(contr.names[i], temp, gs.annot, 
                    file.name[i])
            generateAllHeatmapsPage(contr.names[i], temp, gs.annot, 
                    file.name.hm[i])
            if (length(grep("^kegg", gs.annot@label)) == 1){
                generateAllPathwaysPage(contr.names[i], temp, 
                        gs.annot, file.name.pv[i])
            }
            if ((gs.annot@label == "c5" || gs.annot@label == "gsdbgo") 
                    && "GOID" %in% colnames(gs.annot@anno) ){
                generateAllGOgraphsPage(contr.names[i], gs.annot, sort.by,
                        file.name.go[i])
            }
            
            generateSummaryPage(contr.names[i], gs.annot, sum.plot.axis, 
                    sort.by, length(egsea.results), file.name.sum[i])
        }   
        # print out the comparison analysis pages
        if (length(gsets.top) > 1){      
            egsea.comparison.all = egseaResults@results[[gs.annot@label]]$comparison[["test.results"]]
            egsea.comparison.all = egsea.comparison.all[order(egsea.comparison.all[,sort.by]), ]
            egsea.comparison = egsea.comparison.all[1 : ifelse(
                                nrow(egsea.comparison.all) > display.top, 
                                display.top, nrow(egsea.comparison.all)), ]           
            # Comparison analysis reports generated here
            generateSumPlots.comparison(egsea.results = egsea.results, 
                    egsea.comparison = egsea.comparison.all, 
                    gs.annot = gs.annot, gsa.dir = egsea.dir, 
                    sum.plot.cutoff=sum.plot.cutoff, 
                    sum.plot.axis=sum.plot.axis)    
            plotHeatMapsLogFC.comparison(gene.sets = 
                            rownames(egsea.comparison), 
                    fc=logFC,
                    limma.tops = limma.tops,
                    gs.annot=gs.annot, symbolsMap=symbolsMap,
                    gsa.dir=egsea.dir)
            file.name.bar = paste0(summary.dir,                      
                    "comparison-", gs.annot@label, "-bar-plot-",
                    sort.by) 
            if (!file.exists(paste0(file.name.bar, ".png")))
                plotBars(egseaResults, gs.annot@label, contrast = 0, 
                    number = display.top, sort.by = sort.by,
                    bar.vals = "p.adj",
                    file.name = file.name.bar, 
                    format = NULL, verbose=verbose)
            
            file.name = paste0(ranked.gs.dir, "/ranked-", gs.annot@label, 
                    "-gene-sets-compare.txt")   
            file.name.hm = paste0(egsea.dir, "/hm-top-gs-", gs.annot@label, 
                    "/allHeatmaps.html")                
            file.name.sum = paste0(summary.dir,gs.annot@label, 
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
                generateAllGOgraphsPage.comparison(contr.names, 
                        gs.annot, sort.by, file.name.go) 
            }
            
            generateSummaryPage.comparison(contr.names, gs.annot, sum.plot.axis,
                    sort.by, file.name.sum)  
        }         
    }    
    createEGSEAReport(egseaResults@sampleSize, contr.names, 
            gs.annots, egseaResults@baseInfo, 
            combineMethod, sort.by,  egsea.dir,
            logFC.calculated, symbolsMap,
            egseaResults@egsea.version,
            egseaResults@egseaData.version)    
    if (interactive()) try(browseURL(paste0("file://", 
                                normalizePath(egsea.dir),
                                '/index.html')))
    elapsed.time = proc.time() - start.time    
    timestamp()
    print(paste0("EGSEA report generation took ", elapsed.time["elapsed"], " seconds."))
    print("EGSEA report has been generated.")
}

write.base.stats <- function(base.results, baseGSEAs, contr.names,
        gs.annot, ranked.gs.dir, num.threads){
    args.all = list()
    i = 1
    for (method in baseGSEAs){
        for (contrast in contr.names){                   
            base.results.ind = base.results[[contrast]][[method]]
            file.name = paste0(ranked.gs.dir, "/", method,
                    "-ranked-", 
                    gs.annot@label, "-gene-sets-", 
                    contrast, '.txt')                    
            args.all[[i]] = list(
                    base.frame = base.results.ind[order(
                                    base.results.ind[,"Rank"]), ],
                    contrast = contrast,
                    gs.annot = gs.annot,
                    method = toupper(method),
                    file.name = file.name)
            i = i + 1
        }
    }
    if (Sys.info()['sysname'] == "Windows" ||
            length(args.all) <= 1)
        results = lapply(args.all, write.base.results)
    else
        results = mclapply(args.all, write.base.results, 
                mc.cores = num.threads)
    
}
write.base.results <- function(args){
    writeResultsToHTML(args$contrast, args$base.frame, 
            args$gs.annot, args$method, args$file.name)
}