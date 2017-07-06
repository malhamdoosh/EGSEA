# Wrapper function to run PADOG on different contrasts

runpadog <- function(voom.results, contrast, gs.annot, 
        num.workers=4, verbose = TRUE){     

    # run padog and write out ranked 'gene sets' for each 'contrast'    
    if (is.matrix(contrast)){
        contr.names = colnames(contrast)
        contr.num = ncol(contrast)
    }else{
        contr.names = names(contrast)
        contr.num = length(contrast)
    }
    
    groupData = prepareTwoGroupsData(voom.results, contrast, gs.annot,
            min.samples = 3, verbose)
#    padog.results = vector("list", ncol(contrast))
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
        padog.results = lapply(args.all, runpadog.contrast)
    else
        padog.results = mclapply(args.all, runpadog.contrast, 
                mc.cores=num.workers)
#    names(padog.results) = colnames(contrast)
    return(padog.results)
}

runpadog.contrast <- function(args){
    if (args$verbose)
        message("   Running PADOG for ", args$contrast)
    else
        message(".", appendLF = FALSE)
    
    group = c(rep("c", length(args$group1)), 
            rep("d", length(args$group2)))
    
    padog.results = padog(esetm=args$logCPM, 
            group=group, paired=FALSE, 
            gslist=args$gsets, NI=100, 
            verbose=FALSE)
    
    padog.results = padog.results[order(
                    padog.results[,"Ppadog"], 
                    -padog.results[,"padog0"]),]   
    padog.results = cbind(
            Rank = seq(1, nrow(padog.results)), 
            padog.results)
    colnames(padog.results)[which(colnames(padog.results) == "Ppadog")] = "p.value"  
    padog.results = cbind(padog.results,
            p.adj=p.adjust(padog.results[, "p.value"], method="BH"))
    return(padog.results)
}



