# Wrapper function to run GAGE on different contrasts

rungage <- function(voom.results, contrast, gs.annot, 
        num.workers=4, verbose = TRUE){     
    # run GAGE and write out ranked 'gene sets' for each 'contrast'
    if (is.matrix(contrast)){
        contr.names = colnames(contrast)
        contr.num = ncol(contrast)
    }else{
        contr.names = names(contrast)
        contr.num = length(contrast)
    }
      
    groupData = prepareTwoGroupsData(voom.results, contrast, gs.annot,
            min.samples = 3, verbose)
    #gage.results = vector("list", ncol(contrast))    
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
        gage.results = lapply(args.all, rungage.contrast)
    else
        gage.results = mclapply(args.all, rungage.contrast, 
                mc.cores=num.workers)
    #stop("here")
    #names(gage.results) = colnames(contrast)
    return(gage.results)
}

rungage.contrast <- function(args){
    if (args$verbose)
        message("   Running GAGE for ", args$contrast)
    else
        message(".", appendLF = FALSE)
    groupData = args$groupData
    # same.dir=FALSE ==> Two directional test
    gage.results = gage(exprs=args$logCPM, 
            gsets=args$gsets, 
            ref=args$group1,
            samp=args$group2, same.dir=FALSE, 
            compare = "unpaired")$greater[, 1:5]        
    # returns PropDown/PropUp ==> proportion of genes that are 
# down/up-regulated
    gage.results = gage.results[ order ( 
                    gage.results[,"p.val"]),]  
    gage.results = cbind(Rank=seq(1, nrow(gage.results)), 
            gage.results)    
    colnames(gage.results)[which(colnames(gage.results) == "p.val")] = "p.value"
    colnames(gage.results)[which(colnames(gage.results) == "q.val")] = "p.adj"
    return(gage.results)
}


