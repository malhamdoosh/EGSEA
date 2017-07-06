# Wrapper function to run GLOBALTEST on different contrasts


runglobaltest <- function(voom.results, contrast, gs.annot, 
        num.workers=4, verbose = TRUE){     

    # run globaltest and write out ranked 'gene sets' for each 'contrast'  
    if (is.matrix(contrast)){
        contr.names = colnames(contrast)
        contr.num = ncol(contrast)
    }else{
        contr.names = names(contrast)
        contr.num = length(contrast)
    }
    
    gt.options(transpose=TRUE)     
    groupData = prepareTwoGroupsData(voom.results, contrast, gs.annot, 
            verbose = verbose)    
#    globaltest.results = vector("list", ncol(contrast))
    args.all = list()
    for(i in 1:contr.num){
        args.all[[i]] = list(contrast = contr.names[i],
                logCPM = groupData$data[[i]]$logCPM,
                group1 = groupData$data[[i]]$group1,
                group2 = groupData$data[[i]]$group2,
                gsets = groupData$gsets,
                gs.annot = gs.annot,                
                verbose = verbose)
    }
    names(args.all) = contr.names
    if (Sys.info()['sysname'] == "Windows" || contr.num <= 1)
        globaltest.results = lapply(args.all, runglobaltest.contrast)
    else
        globaltest.results = mclapply(args.all, runglobaltest.contrast, 
                mc.cores=num.workers)
#    names(globaltest.results) = colnames(contrast)
    return(globaltest.results)
}

runglobaltest.contrast <- function(args){
    if (args$verbose)
        message("   Running GLOBALTEST for ", 
                        args$contrast)
    else
        message(".", appendLF = FALSE)
    
    groupIndx = c(args$group2, 
            args$group1)
    data.log.sel = args$logCPM[, groupIndx]
    group = c(rep("d", length(args$group2)), 
            rep("c", length(args$group1)))
    colnames(data.log.sel) = group
    # perform the globaltest
    gtobject = gt(factor(group), 
            data.log.sel, subsets=args$gs.annot@idx, 
            permutations=10000)
#    print("heregt")
    globaltest.results = result(gtobject)
#    print(head( globaltest.results))
    globaltest.results = globaltest.results[
            order(globaltest.results[,"p-value"],
                    -globaltest.results[,"Statistic"]),]
    globaltest.results = cbind(
            Rank=seq(1, nrow(globaltest.results)),             
            globaltest.results)
    colnames(globaltest.results)[which(
                    colnames(globaltest.results) == "p-value")] = "p.value"
    padog.results = cbind(globaltest.results,
            p.adj=p.adjust(globaltest.results[, "p.value"], method="BH"))
#    print(head( globaltest.results))
    return(globaltest.results)
}

