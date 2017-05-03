# Wrapper function to run FRY on different contrasts

runfry <- function(voom.results, contrast, gs.annot,
        num.workers=4, verbose = TRUE){     
    # run fry and write out ranked 'gene sets' for each 'contrast'    
    if (is.matrix(contrast)){
        contr.names = colnames(contrast)
        contr.num = ncol(contrast)
    }else{
        contr.names = names(contrast)
        contr.num = length(contrast)
    }
     
#    fry.results = vector("list", ncol(contrast))  
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
        fry.results = lapply(args.all, runfry.contrast)
    else
        fry.results = mclapply(args.all, runfry.contrast, 
                mc.cores=num.workers)    
#    names(fry.results) = colnames(contrast)
    return(fry.results)
}

runfry.contrast <- function(args){
    if (args$verbose)
        print(paste0("   Running FRY for ", args$contrast.name))
    else
        cat(".")
    capture.output(fry.results <- fry(y=args$voom.results, 
                    index=args$gs.annot@idx, 
                    design=args$voom.results$design, 
                    contrast=args$contrast))
    # returns PropDown/PropUp ==> proportion of genes that are 
# down/up-regulated
    fry.results = fry.results[order(fry.results[, 
                            "PValue"]),]
    fry.results = cbind(Rank=seq(1, 
                    nrow(fry.results)), fry.results)
    colnames(fry.results)[which(colnames(fry.results) == "PValue")] = "p.value"
    colnames(fry.results)[which(colnames(fry.results) == "FDR")] = "p.adj"
    return(fry.results)
}

#y <- matrix(rnorm(100*4),100,4)
#design <- cbind(Intercept=1,Group=c(0,0,1,1))
## First set of 5 genes contains 3 that are genuinely differentially expressed
#index1 <- 1:5
#y[index1,3:4] <- y[index1,3:4]+3
## Second set of 5 genes contains none that are DE
#index2 <- 6:10
#
#r <- fry(y,list(set1=index1,set2=index2),design,contrast=2)