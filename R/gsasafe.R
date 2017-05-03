# Wrapper function to run SAFE (Significance Analysis of Function and 
# Expression) on different contrasts

runsafe <- function(voom.results, contrast, gs.annot, 
        num.workers=4, verbose = TRUE){     

    # run safe and write out ranked 'gene sets' for each 'contrast' 
    if (is.matrix(contrast)){
        contr.names = colnames(contrast)
        contr.num = ncol(contrast)
    }else{
        contr.names = names(contrast)
        contr.num = length(contrast)
    }
    
    groupData = prepareTwoGroupsData(voom.results, contrast, gs.annot,
            min.samples = 3, verbose)
    capture.output(C.mat <- getCmatrix(keyword.list=groupData$gsets, 
                    present.genes=rownames(groupData$data[[1]]$logCPM)))
#    safe.results = vector("list", ncol(contrast))  
    args.all = list()
    for(i in 1:contr.num){
        args.all[[i]] = list(contrast = contr.names[i],
                C.mat = C.mat,
                logCPM = groupData$data[[i]]$logCPM,
                group1 = groupData$data[[i]]$group1,
                group2 = groupData$data[[i]]$group2,
                gsets = groupData$gsets,
                gs.annot = gs.annot,                
                verbose = verbose)
    }
    names(args.all) = contr.names
    if (Sys.info()['sysname'] == "Windows" || contr.num <= 1)
        safe.results = lapply(args.all, runsafe.contrast)
    else
        safe.results = mclapply(args.all, runsafe.contrast, 
                mc.cores=num.workers)
#    names(safe.results) = colnames(contrast)
    return(safe.results)
}

runsafe.contrast <- function(args){
    if (args$verbose)
        print(paste0("   Running SAFE for ", args$contrast))
    else
        cat(".")
    group = c(rep("Ctr", length(args$group1)), 
            rep("Trt", length(args$group2)))
    
    capture.output(temp <- safe(X.mat=args$logCPM, 
                    y.vec=group, C.mat=args$C.mat, 
                    print.it=FALSE, Pi.mat=100))#, 
    safe.results = safe.toptable(temp, 
            number=length(args$gsets),
            description = FALSE)
    safe.results[, "P.value"] = as.numeric(
            safe.results[, "P.value"])
    safe.results[, "Adj.p.value"] = as.numeric(
            safe.results[, "Adj.p.value"])
    rownames(safe.results) = safe.results[, "GenesetID"]
    safe.results = safe.results[, -c(1,6)]        
    safe.results = safe.results[order(
                    safe.results[,"P.value"],
                    -safe.results[,"Statistic"]),]
    safe.results = cbind(Rank = seq(1, 
                    nrow(safe.results)), safe.results)
    colnames(safe.results)[which(colnames(safe.results) == "P.value")] = "p.value"
    colnames(safe.results)[which(colnames(safe.results) == "Adj.p.value")] = "p.adj"
    return(safe.results)
}

