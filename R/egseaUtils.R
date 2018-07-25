#Ensemble of Gene Set Enrichment Analyses 
# 
# Author: Monther Alhamdoosh, E:m.hamdoosh@gmail.com


egsea.main <- function(voom.results, contrast, gs.annots, baseGSEAs, 
                combineMethod, combineWeights, sort.by, report.dir, 
                kegg.dir, logFC, symbolsMap, minSize, display.top, 
                logFC.cutoff, fdr.cutoff, sum.plot.cutoff, sum.plot.axis, 
                vote.bin.width, keep.base, verbose, num.threads, 
                report, keep.limma, keep.set.scores, interactive){    
    message("EGSEA analysis has started")
    start.time <- proc.time()
    timestamp()
    # check arguments are valid    
#    print(class(voom.results))
    stopifnot((class(voom.results) == "list" && 
                        "ids" %in% names(voom.results)) 
                    || class(voom.results) == "EList")  
    #stopifnot(!is.null(contrast))
    stopifnot(!is.null(gs.annots))
    stopifnot(length(baseGSEAs) > 0 && length(setdiff(baseGSEAs, egsea.base())) == 0)
    stopifnot(combineMethod %in% egsea.combine())
    stopifnot(sort.by %in% c(egsea.sort()[1:8], baseGSEAs))
    baseGSEAs = unique(sapply(baseGSEAs, tolower))    
    combineMethod = tolower(combineMethod)  
    sort.by = tolower(sort.by)   
    if (length(baseGSEAs) <= 1){
        warning("The ensemble mode was disabled. No sufficient base methods.")
        if (! sort.by %in% c("p.adj", "p.value")){
            sort.by = "p.adj"
            warning("The argument 'sort.by' was set to \"p.adj\"")
        }
    }
    keep.set.scores = (keep.set.scores && length(intersect(baseGSEAs, 
                    c("ssgsea")))) # , "gsva", "plage", "zscore"  
    if (is.null(sum.plot.axis))
        sum.plot.axis = "p.adj"
    # create a list of GSCollectionIndex objects
    if (class(gs.annots) == "GSCollectionIndex"){
        gs.annot = gs.annots
        gs.annots = list()
        gs.annots[[gs.annot@label]] = gs.annot
    }  
    # format the contrast matrix/vector and create contrast names
    if (is.null(contrast)){
        if (is.null(voom.results$targets) || 
                is.null(voom.results$targets$group))
            stop(paste0("The data frame 'targets' of the object 'voom.results' ",
                            "must have a column named 'group'."))
        group.levels = levels(factor(voom.results$targets$group))
        if (all(group.levels %in% colnames(voom.results$design))){
            contrast = makeContrastPairs(
                    colnames(voom.results$design),
                    group.levels)
        }else
            contrast = 2:length(group.levels)        
    }
    if (!is.matrix(contrast)){
        stopifnot(max(contrast) <= ncol(voom.results$design))
        message("The argument 'contrast' is recommended to be a matrix object.\n",
                    "See Vignette or Help.")
        if (is.null(names(contrast))){
            if (is.null(voom.results$targets) || 
                    is.null(voom.results$targets$group))
                stop("The data frame 'targets' of the object 'voom.results' ",
                                "must have a column named 'group'.")
            ref.group = levels(factor(voom.results$targets$group))[1]  
            names(contrast) = paste0(colnames(voom.results$design)[contrast], 
                    "vs", ref.group)
        }
        contr.names = gsub("[[:punct:]]+", "", names(contrast))
        contr.names = gsub("[[:space:]]+", "", names(contrast))
        names(contrast) = contr.names
        colnames(voom.results$design)[contrast] = contr.names
    }else {
        stopifnot(nrow(contrast) == ncol(voom.results$design))
        if (is.null(colnames(contrast)))
            colnames(contrast) = paste0("Contrast", rep(1, ncol(contrast)))
        contr.names = gsub("[[:punct:]]", "", colnames(contrast))
        contr.names = gsub("[[:space:]]+", "", colnames(contrast))
        colnames(contrast) = contr.names
    }   
    # calculate logFC from limma DE analysis 
    logFC.calculated = "No"
    if (is.null(logFC)){
        #row names should be Entrez Gene IDs in order to plot KEGG pathways
        message("Log fold changes are estimated using limma package ... ")
        tmp = runStandardLimmaDEA(voom.results, contrast, logFC.cutoff, fdr.cutoff)
        logFC = tmp$logFC
        limma.results = tmp$limma.results
        limma.tops = tmp$limma.tops
        logFC.calculated = "Yes"
    }else if (!is.matrix(logFC) || !identical(colnames(logFC), contr.names)){
        stop(paste0("logFC should be a matrix object with column names equal ", 
                        "to the (column) names of the argument 'contrast'."))
    }else if (class(voom.results) == "EList"){
        tmp = runStandardLimmaDEA(voom.results, contrast, logFC.cutoff, fdr.cutoff)
        limma.results = tmp$limma.results
        limma.tops = tmp$limma.tops
    }else{
        limma.results = new("MArrayLM")
        limma.tops = list()
    }    
    # check featureIDs and logFC matrix row names
    featureIDs = as.character(gs.annots[[1]]@featureIDs)
    if (!identical(rownames(logFC), featureIDs)){     
        logFC = logFC[match(featureIDs, rownames(logFC)) , ]    
        if (!identical(rownames(logFC), featureIDs)){
            stop("The row names of the fold change matrix should \n", 
                  "match the featureIDs vector in the gs.annot list.\n",
                 paste(rownames(logFC)[1:10], collapse=' '), "\n", 
                 paste(featureIDs[1:10], collapse= ' '))          
        }
    }    
    # check the 'symbolsMap' argument and replace NA symbols 
    if (is.null(symbolsMap)){
        symbolsMap = data.frame()
    }else if (nrow(symbolsMap) > 0 && ncol(symbolsMap) >= 2){
        na.sym = is.na(symbolsMap[, 2])
        if (sum(na.sym) > 0){
            warning("Some \"NA\" Gene Symbols were replaced with Feature IDs")
            symbolsMap[na.sym, 2] = symbolsMap[na.sym, 1]
        }
    }    
    if (is.null(report.dir))
      report.dir = ""
    # optimize number of cores to be used 
    num.threads = optimizeNumThreads(num.threads, length(baseGSEAs), 
            length(contr.names), verbose)    
    # create the EGSEAResults object to be populated
    gsas = EGSEAResults(contr.names = contr.names, 
            contrast = contrast,
            sampleSize = getNumberofSamples(voom.results, contrast), 
            gs.annots = gs.annots, baseMethods=baseGSEAs, 
            baseInfo = getBaseInfo(baseGSEAs),
            combineMethod = combineMethod, sort.by = sort.by,   
            symbolsMap = symbolsMap,
            logFC = logFC, logFC.calculated = logFC.calculated,
            sum.plot.axis = sum.plot.axis, 
            sum.plot.cutoff = sum.plot.cutoff,
            report = report, report.dir = report.dir,
            egsea.version = paste0(packageVersion("EGSEA")),
            egseaData.version = paste0(packageVersion("EGSEAdata"))
            )  
    if (keep.limma && length(limma.results) > 0)
        gsas@limmaResults = limma.results
    ###### START THE EGSEA ANALYSIS ON INDIVIDUAL COLLECTIONS ######
    skipped = c()    
    for (gs.annot in gs.annots){
        gs.annot = selectGeneSets(gs.annot, min.size=minSize)  
        if (length(gs.annot@idx) == 0){
            message("No gene sets in ", gs.annot@label, 
                    " meets the minimum size criterion.")
            skipped = c(skipped, gs.annot@label)
            next
        }
        # run egsea and write out ranked gene sets for all contrasts
        message("EGSEA is running on the provided data and ",
                        gs.annot@label, " collection")   
        
        results <- runegsea(voom.results = voom.results, 
                contrast = contrast, limma.tops = limma.tops,
                baseGSEAs = baseGSEAs,
                combineMethod = combineMethod, combineWeights = NULL,
                gs.annot = gs.annot, logFC = logFC, 
                logFC.cutoff = logFC.cutoff, fdr.cutoff = fdr.cutoff,
                vote.bin.width=vote.bin.width, keep.base = keep.base,
                num.workers = num.threads, verbose=verbose)            
             
        egsea.results = results[["egsea.results"]]
        # order results based on the sort.by argument
        for (i in 1:length(egsea.results)){
            # sort based on the average ranking
            egsea.results[[i]] = egsea.results[[i]][
                                order(egsea.results[[i]][,sort.by],
                                decreasing=(sort.by == "significance")), 
                                ]      
        }  
        
        gsa = list("test.results"=egsea.results)
        if (keep.base)
            gsa[["base.results"]] = results[["base.results"]]        
        if (keep.set.scores){
            set.methods = intersect(baseGSEAs, 
                    c("ssgsea")) # , "gsva", "plage", "zscore"
            gsa[["set.scores"]] = calculateSetScores(voom.results, gs.annot, 
                    set.methods, num.threads, verbose)
        }  
        # Comparison analysis reports generated here
        if (length(egsea.results) > 1){         
            egsea.comparison = createComparison(egsea.results, 
                    combineMethod = combineMethod, 
                    display.top=Inf, sort.by = sort.by)
            gsa[["comparison"]] = list()
            gsa$comparison[["test.results"]] = egsea.comparison
            egsea.comparison.all = egsea.comparison            
            egsea.comparison = egsea.comparison[1:ifelse(nrow(egsea.comparison) 
                > display.top, 
                display.top, nrow(egsea.comparison)), ]
            # gsa$comparison[["top.gene.sets"]] = rownames(egsea.comparison)
        }   
        gsas = addEGSEAResult(gsas, gs.annot@label, gsa)
        #gsas[[gs.annot@label]] = gsa
    }    
    elapsed.time = proc.time() - start.time    
    timestamp()
    message("EGSEA analysis took ", elapsed.time["elapsed"], " seconds.")
    message("EGSEA analysis has completed")
    if (report){
      generateMainReport(gsas, limma.tops = limma.tops,
               display.top = display.top,
               kegg.dir = kegg.dir, num.threads = num.threads, 
               print.base = TRUE,
               interactive = interactive,
               verbose = verbose) 
    }
    return(gsas)
}

optimizeNumThreads <- function(num.threads, base.num, contr.num, verbose=TRUE){    
    cores = detectCores()    
    exp.threads.num = min(num.threads, base.num) * 
                    min(num.threads, contr.num)
    if (verbose)
        message("Expected number of running processes: ", 
                        exp.threads.num)
    if (cores < exp.threads.num ){
        new.threads.num = 1
        exp.threads.num = min(new.threads.num, base.num) * 
                min(new.threads.num, contr.num)
        while (new.threads.num < cores &&
                exp.threads.num  < cores ){
            new.threads.num = new.threads.num + 1
            exp.threads.num = min(new.threads.num, base.num) * 
                    min(new.threads.num, contr.num)
#            print(new.threads.num)
#            print(exp.threads.num)
        }        
    } else if (cores < num.threads){
        new.threads.num = cores
    }
    else
        new.threads.num = num.threads
    if (new.threads.num  != num.threads){      
        message("Number of used cores has changed to ", 
                        new.threads.num , "\nin order to avoid ",
                        "CPU overloading.")
    }
    return(new.threads.num)
}

calculateSetScores <- function(voom.results, gs.annot, set.methods,
        num.workers, verbose = TRUE){    
    set.scores = list()
    if (length(set.methods) == 0){
        message("The parameter keep.set.scores has nothing to keep.")
        return(set.scores)
    }
    data.log = voom.results$E
    rownames(data.log) = as.character(seq(1, nrow(data.log)))
    gsets = list()
    for (j in 1:length(gs.annot@idx)){
        gsets[[j]] = as.character(gs.annot@idx[[j]])
    }
    names(gsets) = names(gs.annot@idx)
    for (method in set.methods){
        if (verbose)
            message("Gene set enrichment scores per sample are\n", 
                    "being calculated using ", 
                        method, "...")
        gs.es = calculateSetScores.parallel(data.log, gsets, method, num.workers)
        set.scores[[method]] = gs.es
    }
    return(set.scores)
}

getBaseInfo <- function(baseGSEAs){
    basePkg = list()
    basePkg[["camera"]] = "limma"
    basePkg[["roast"]] = "limma"
    basePkg[["safe"]] = "safe"
    basePkg[["gage"]] = "gage"
    basePkg[["padog"]] = "PADOG"
    basePkg[["plage"]] = "GSVA"
    basePkg[["zscore"]] = "GSVA"
    basePkg[["gsva"]] = "GSVA"
    basePkg[["ssgsea"]] = "GSVA"
    basePkg[["globaltest"]] = "globaltest"
    basePkg[["ora"]] = "stats"
    basePkg[["fry"]] = "limma"    
    baseInfo = list()
    for (baseGSEA in baseGSEAs){
        baseInfoi = list()
        pkg = basePkg[[baseGSEA]]
        baseInfoi$version = packageVersion(pkg)        
        baseInfoi$package = pkg
        baseInfo[[baseGSEA]] = baseInfoi
    }
    return(baseInfo)
}


egsea.selectTopGeneSets <- function(egsea.results, display.top, gs.annot, 
 verbose=TRUE){
    contrast.names = names(egsea.results)
    top.gene.sets = list()
    
    for(i in 1:length(egsea.results)){      
        #num.gene.sets.fdr = sum(egsea.results[[i]]$FDR < fdr[i])
        num.gene.sets.fdr = ifelse(length(gs.annot@idx) > display.top, 
display.top,length(gs.annot@idx))
        top.print = ifelse(length(gs.annot@idx) >= 10, 10,length(gs.annot@idx))
        if (top.print > num.gene.sets.fdr)
            top.print = num.gene.sets.fdr
        if (num.gene.sets.fdr > 0){
            egsea.results.top = egsea.results[[i]][1:num.gene.sets.fdr,]    
        
            top.gene.sets[[i]] = rownames(egsea.results.top) 
            if (verbose){
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
                                        gs.annot@anno[,2]), -2],
                            egsea.results.top)
                    print(top.table[1:top.print, c("Type", "p.adj")])
                }                
            }
            
        }else{            
            message("No gene sets found")
            top.gene.sets[[i]] = NA
        }
        
    }
    
    names(top.gene.sets) = contrast.names
    return(top.gene.sets)
}


runbaseGSEAParallelWorker <- function(args){
    #print(paste0("Running ", toupper(args$baseGSEA), " on all contrasts ... "))
    tryCatch({
                #t = system.time(
                temp.result <- runbaseGSEA(method=args$baseGSEA, 
                        args$voom.results, args$contrast, args$gs.annot,                        
                        num.threads = args$num.threads, verbose=args$verbose)
                #)
                #print(paste0(args$baseGSEA, " ", args$gs.annot$name, " ", t[3]))
                if (args$verbose)
                    message(paste0("Running ", toupper(args$baseGSEA), " on all ", 
							"contrasts ... COMPLETED"))
                else{
                    message(args$baseGSEA, "*", appendLF = FALSE)
                }
                return(temp.result)
            }, 
            error = function(e) {
                message(toupper(args$baseGSEA), 
                           " encountered an error -> ", e )
            })
    return(NULL)
}

runegsea <- function(voom.results, contrast, limma.tops,
            baseGSEAs, combineMethod, 
            combineWeights=NULL, gs.annot, logFC, 
            logFC.cutoff, fdr.cutoff, 
            vote.bin.width, keep.base=TRUE, 
            num.workers=8, verbose=FALSE){     
    stopifnot(class(gs.annot) == "GSCollectionIndex")
    if (!is.null(voom.results$E))
        geneIDs = rownames(voom.results$E)
    else 
        geneIDs = voom.results$ids
    stopifnot(identical(geneIDs, gs.annot@featureIDs))
    if (is.matrix(contrast)){
        contr.names = colnames(contrast)        
    }else{
        contr.names = names(contrast)        
    }
    # egsea.results.details stores Contrasts ==> Individual Results
    egsea.results.details = vector("list", length(contr.names))  
    names(egsea.results.details) = contr.names
    for (i in 1:length(contr.names)){
        egsea.results.details[[i]] = vector("list", length(baseGSEAs))
        names(egsea.results.details[[i]]) = baseGSEAs
    }
    # run enrichment analysis using base methods
    args.all = list()
    threads.per.base = num.workers / 2 # ceiling(num.workers / length(baseGSEAs))
    for (baseGSEA in baseGSEAs){
        args.all[[baseGSEA]] = list(baseGSEA=baseGSEA, 
                    voom.results=voom.results,
                    contrast=contrast, gs.annot=gs.annot,                                        
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
        message("")
    # collect results   
    #print(names(temp.results))
    for (baseGSEA in baseGSEAs){
        for (i in 1:length(contr.names)){
            # order is important when combine
            if (is.null(temp.results[[baseGSEA]])){
                stop("ERROR: One of the base methods failed on this ", 
                "dataset (", baseGSEA , 
                ").\nRemove it and try again.\nSee error messages for ", 
                "more information.")
            }
#            print(paste0(baseGSEA, colnames(contrast)[i]))
#            if (baseGSEA == "globaltest")
#                print(temp.results[[baseGSEA]][[i]])
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
        gs.avg.fcs.dir = numeric(0)
        gs.dirs = numeric(0)
        gsets = as.character(rownames(egsea.results[[i]]))
        fc = logFC[, i]
        if (length(limma.tops) > 0)
            t = limma.tops[[i]]        
        for (j in 1:length(gsets)){             
            sel.genes = gs.annot@idx[[gsets[j]]]
            gset.fc = fc[sel.genes]    
            if (length(limma.tops) > 0){
                temp = gset.fc[t[sel.genes, "adj.P.Val"] <= fdr.cutoff]
                if (length(temp) > 0)
                    gset.fc = temp            
                temp = abs(gset.fc)
                temp = temp[temp >= logFC.cutoff]
                if (length(temp) == 0)
                    temp = abs(gset.fc)
            }else
                temp = abs(gset.fc)
            gs.avg.fcs = c(gs.avg.fcs, mean(temp, na.rm=TRUE))          
            up = sum(gset.fc > logFC.cutoff, na.rm=TRUE)
            dn = sum(gset.fc < -logFC.cutoff, na.rm=TRUE)
            if (up > dn){
                gs.dirs = c(gs.dirs, 1)
                gs.avg.fcs.dir = c(gs.avg.fcs.dir, 
                        mean(gset.fc[gset.fc > logFC.cutoff], na.rm=TRUE))
            }else{
                gs.dirs = c(gs.dirs, -1)
                gs.avg.fcs.dir = c(gs.avg.fcs.dir, 
                        mean(gset.fc[gset.fc < -logFC.cutoff], na.rm=TRUE))
            }            
        }   
#        print(head(egsea.results[[i]]))
        pvalues = egsea.results[[i]][, "p.adj"]        
        pvalues = -1 * log10(pvalues)
        pvalues[pvalues == Inf] = max(pvalues[pvalues != Inf]) + 50
        pvalues[is.na(pvalues)] = 0 
        sig = pvalues * gs.avg.fcs
        if (max(sig, na.rm=TRUE) != min(sig, na.rm=TRUE))
            sig = (sig - min(sig, na.rm=TRUE)) / (max(sig, na.rm=TRUE) - 
                    min(sig, na.rm=TRUE)) * 100
        m = length(baseGSEAs)       
        if (m == 1){            
            egsea.results[[i]] = cbind(egsea.results[[i]], 
                        "avg.logfc" = gs.avg.fcs,
                        "avg.logfc.dir" = gs.avg.fcs.dir,
                    "direction" = gs.dirs, "significance" = sig)
        }
        else{
            n = ncol(egsea.results[[i]])
            # insert stat columns in the middle
            egsea.results[[i]] = cbind(egsea.results[[i]][, 1:(n-m)], 
                            "avg.logfc" = gs.avg.fcs,
                            "avg.logfc.dir" = gs.avg.fcs.dir,
                    "direction" = gs.dirs, "significance" = sig, 
                    egsea.results[[i]][, (n-m+1):n] )
        }        
    }
    
    names(egsea.results) = contr.names
    results = list("egsea.results"=egsea.results)
    if (keep.base){
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
        } else if (combineMethod == "votep"){        
          data[data == 0] = 1*10^-22
          data[data == 1] = 1 - 1*10^-5
          pvalues = sapply(apply(data,  1, function(y) votep(y[!is.na(y)])), 
                           function(x) x$p)            
        } else if (combineMethod == "median"){        
          data[data == 0] = 1*10^-22
          data[data == 1] = 1 - 1*10^-5
          pvalues = sapply(apply(data,  1, function(y) median(y[!is.na(y)])), 
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
    if (length(results.multi) == 1 && 
        length(results.multi[[1]]) == 1 && 
        names(results.multi[[1]]) == "ora"){
        results.multi[[1]] = results.multi[[1]][[1]]
        results.multi[[1]] = results.multi[[1]][, 
                colnames(results.multi[[1]]) != "Rank"]
        return(results.multi)
    }
  
    if (length(results.multi[[1]]) == 1){
        results.combined = vector('list', length(results.multi))
        names(results.combined) = names(results.multi)
        for (i in 1:length(results.combined)){
            results.combined[[i]] = results.multi[[i]][[1]]
            results.combined[[i]] = results.combined[[i]][, 
                            colnames(results.combined[[i]]) != "Rank"]
        }
        return(results.combined)
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
            results.comp[[i]][,j] = as.numeric(results.multi[[i]][[j]][, 
                        "p.adj"])           
            results.ranked[[i]][,j] = as.numeric(results.multi[[i]][[j]][, 
                            "Rank"])
            if (max(results.comp[[i]][,j], na.rm = TRUE) < 0.000001)
                warning(method, 
                    " produces very low p-values on ", 
                    names(results.comp)[i])
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
            results.comp[[i]][,j] = as.numeric(results.multi[[i]][[j]][, 
                            "p.value"]) 
            results.ranked[[i]][,j] = as.numeric(results.multi[[i]][[j]][, 
                "Rank"])
            if (max(results.comp[[i]][,j], na.rm = TRUE) < 0.000001)
                warning(method, 
                      " produces very low p-values on ", 
                      names(results.comp)[i])
            j = j + 1
        }
    }
    return(list(pvalues=results.comp, ranks=results.ranked))
}

runbaseGSEA <- function(method, voom.results, contrast, gs.annot,                         
                        num.threads = 4, verbose=TRUE){
    if (method == "camera"){        
        return(runcamera(voom.results = voom.results, contrast = contrast, 
            gs.annot = gs.annot, num.workers = num.threads, verbose = verbose))     
    }else if (method == "roast"){       
        return(runroast(voom.results = voom.results, contrast = contrast, 
            gs.annot = gs.annot, num.workers = num.threads, verbose = verbose))     
    }else if (method == "fry"){
        return(runfry(voom.results = voom.results, contrast = contrast, 
            gs.annot = gs.annot, num.workers = num.threads, verbose = verbose))
    }else if (method == "gage"){        
        return(rungage(voom.results = voom.results, contrast = contrast, 
            gs.annot = gs.annot, num.workers = num.threads, verbose = verbose))     
    }else if (method == "padog"){       
        return(runpadog(voom.results = voom.results, contrast = contrast, 
            gs.annot = gs.annot, num.workers = num.threads, verbose = verbose))
    }else if (method == "plage"){       
        return(rungsva(method="plage", voom.results = voom.results, contrast = 
            contrast, gs.annot = gs.annot, num.workers = num.threads, 
            verbose = verbose))
    }else if (method == "zscore"){  
        return(rungsva(method="zscore", voom.results = voom.results, contrast = 
            contrast, gs.annot = gs.annot, num.workers = num.threads, 
            verbose = verbose))
    }else if (method == "gsva"){        
        return(rungsva(method="gsva", voom.results = voom.results, contrast = 
            contrast, gs.annot = gs.annot, num.workers = num.threads, 
            verbose = verbose))
    }else if (method == "ssgsea"){  
        return(rungsva(method="ssgsea", voom.results = voom.results, contrast = 
            contrast, gs.annot = gs.annot, num.workers = num.threads, 
            verbose = verbose))
    }else if (method == "globaltest"){      
        return(runglobaltest(voom.results = voom.results, contrast = contrast, 
            gs.annot = gs.annot, num.workers = num.threads, verbose = verbose))
    }else if (method == "safe"){        
        return(runsafe(voom.results = voom.results, contrast = contrast, 
            gs.annot = gs.annot, num.workers = num.threads, verbose = verbose))
    }else if (method == "ora"){
        return (runora(voom.results = voom.results, contrast = contrast, 
            gs.annot = gs.annot, num.workers = num.threads, verbose = verbose))         
    }else{
        stop("Method not recognized. Type egsea.base() to see supported base 
methods.")
    }
}

createComparison <- function(egsea.results, combineMethod="fisher", display.top=100, 
        sort.by="p.value"){
    egsea.comparison = numeric(0)
    col.names = colnames(egsea.results[[1]])
    col.names.sel = character(0)
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
            col.names.sel = c(col.names.sel, col.names[i], "p.adj")
        }
        else if (col.names[i] == "med.rank"){
            egsea.comparison = cbind(egsea.comparison, rowMedians(temp, 
                    na.rm=TRUE))
            col.names.sel = c(col.names.sel, col.names[i])
        }
        else if (length(grep("min", col.names[i])) > 0){
            minVals = sapply(1:nrow(temp), function(x) min(temp[x, ], 
                    na.rm=TRUE))
            egsea.comparison = cbind(egsea.comparison, minVals)
            col.names.sel = c(col.names.sel, col.names[i])
        }
        else if (col.names[i] %in% c("avg.rank", "direction", "significance", 
                "avg.logfc", "avg.logfc.dir")){
            egsea.comparison = cbind(egsea.comparison, rowMeans(temp, 
                na.rm=TRUE))
            col.names.sel = c(col.names.sel, col.names[i]) 
        }else if (col.names[i] == "vote.rank"){
            if (length(egsea.results) > 2){
                egsea.comparison = cbind(egsea.comparison, voteRank(temp, 
                    bin.width = -1))
                col.names.sel = c(col.names.sel, col.names[i])
            }else{
                egsea.comparison = cbind(egsea.comparison, rowMeans(temp, 
                            na.rm=TRUE))
                col.names.sel = c(col.names.sel, "vote.rank.mean")
            }
        }else{
            next
        }
    }
    rownames(egsea.comparison) = gset.names
#   print(colnames(egsea.comparison))
#   print(col.names)
    colnames(egsea.comparison) = col.names.sel
    egsea.comparison = egsea.comparison[order(egsea.comparison[,sort.by]), ]
    display.top = ifelse(nrow(egsea.comparison) > display.top, display.top, 
nrow(egsea.comparison))
    return(egsea.comparison[1:display.top, ])
    
}

get.toptables <- function(ebayes.results, contrast){
    if (is.matrix(contrast)){        
        contr.names = colnames(contrast)
        coefs = 1:ncol(contrast)
# between the groups
    }else{
        contr.names = names(contrast)
        coefs = contrast
    }
    limma.tops = list()
    for (i in 1:length(coefs)){
        top.table = topTable(ebayes.results, coef=coefs[i], 
                number=Inf, sort.by="none")        
        rownames(top.table) = rownames(ebayes.results)
        limma.tops[[contr.names[i]]] = top.table
    }
    return(limma.tops)
}

runStandardLimmaDEA <- function(voom.results, contrast, 
        logFC.cutoff, fdr.cutoff){
    # to be changed for gene symbols support
    stopifnot(class(voom.results) == "EList")
    message("limma DE analysis is carried out ... ")
    # fit linear model for each gene using limma package functions
    vfit = lmFit(voom.results, design=voom.results$design) # Fit linear model 
# for each gene given a series of arrays
    if (is.matrix(contrast)){
        vfit = contrasts.fit(vfit, contrast) # make all pair-wise comparisons
        contr.names = colnames(contrast)    
        contr.num =  ncol(contrast)
        coefs = 1:ncol(contrast)
# between the groups
    }else{
        contr.names = names(contrast)    
        contr.num = length(contrast)
        coefs = contrast
    }
    ebayes.results = eBayes(vfit) # compute moderated t-statistics, moderated 
#F-statistic, and log-odds of differential expression by empirical 
#    Bayes moderation of the standard errors towards a common value  
    logFC = matrix(0, nrow(ebayes.results), contr.num)
    limma.tops = list()
    for (i in 1:length(coefs)){
        top.table = topTable(ebayes.results, coef=coefs[i], 
                number=Inf, sort.by="none")
        limma.fc = top.table$logFC      
        names(limma.fc) = rownames(ebayes.results)  
        logFC[, i] = limma.fc
        rownames(top.table) = rownames(ebayes.results)
        limma.tops[[contr.names[i]]] = top.table
        de.genes = top.table[top.table[, "adj.P.Val"] <= fdr.cutoff, ]
        if (nrow(de.genes) == 0)
            warning("It seems the contrast ", 
                    contr.names[i], 
                    " has no DE genes at the selected 'fdr.cutoff'.\n",
                    "The 'fdr.cutoff' was ignored in the calculations.")
        if (nrow(de.genes) > 0){
            de.genes = de.genes[abs(de.genes[, "logFC"]) >= logFC.cutoff, ]
            if (nrow(de.genes) == 0)
                warning("It seems the contrast ", 
                    contr.names[i], 
                    " has no DE genes at the selected 'logFC.cutoff'.\n",
                    "The 'logFC.cutoff' was ignored in the calculations.")
        }
    }
    
    rownames(logFC) = rownames(ebayes.results)
    colnames(logFC) = contr.names
    
    return(list(logFC=logFC, limma.results=ebayes.results, limma.tops=limma.tops))
}

# Adapted from Gordon Smyth https://support.bioconductor.org/p/9228/
makeContrastPairs <- function(
        design.cols,
        group.levels){
    n = length(group.levels)
    stopifnot(identical(group.levels, design.cols[1:n]))
    contr = matrix(0, length(design.cols), choose(n, 2))
    rownames(contr) = design.cols    
    colnames(contr) = 1:choose(n,2)
    k = 0
    for (i in 1:(n-1)){
        for (j in (i+1):n){
            k = k + 1
            contr[j, k] = 1
            contr[i, k] = -1
#            print(contr)
            colnames(contr)[k] = paste0(group.levels[j],
                            "vs", group.levels[i])          
        }
    }        
    return(contr)
}

