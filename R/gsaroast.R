# Wrapper function to run ROAST on different contrasts

runroast <- function(voom.results, contrast, gs.annot, 
        num.workers=4, verbose = TRUE){     
    # run ROAST and write out ranked 'gene sets' for each 'contrast'    
    if (is.matrix(contrast)){
        contr.names = colnames(contrast)
        contr.num = ncol(contrast)
    }else{
        contr.names = names(contrast)
        contr.num = length(contrast)
    }
     
#    roast.results = vector("list", ncol(contrast))  
    args.all = list()
    for(i in 1:contr.num){
        if (is.matrix(contrast))
            args.all[[i]] = list(voom.results = voom.results,
                    contrast.name = contr.names[i],
                    contrast = contrast[,i],
                    gs.annot = gs.annot,
                    verbose = verbose
            )
        else
            args.all[[i]] = list(voom.results = voom.results,
                    contrast.name = contr.names[i],
                    contrast = contrast[i],
                    gs.annot = gs.annot,                    
                    verbose = verbose
            )
    }
    names(args.all) = contr.names
    if (Sys.info()['sysname'] == "Windows" || contr.num <= 1)
        roast.results = lapply(args.all, runroast.contrast)
    else
        roast.results = mclapply(args.all, runroast.contrast, 
                mc.cores=num.workers)
#    names(roast.results) = colnames(contrast)
    return(roast.results)
}

runroast.contrast <- function(args){
    if (args$verbose)
        print(paste0("   Running ROAST for ", args$contrast.name))
    else
        cat(".")
    roast.results = mroast(y=args$voom.results, 
            index=args$gs.annot@idx, 
            design=args$voom.results$design, 
            contrast=args$contrast, nrot=999)
    # returns PropDown/PropUp ==> proportion of genes that are 
# down/up-regulated
    roast.results = roast.results[order(roast.results[, 
                            "PValue"]),]
    roast.results = cbind(Rank=seq(1, 
                    nrow(roast.results)), roast.results)
    colnames(roast.results)[which(colnames(roast.results) == "PValue")] = "p.value"
    colnames(roast.results)[which(colnames(roast.results) == "FDR")] = "p.adj"
    return(roast.results)
}