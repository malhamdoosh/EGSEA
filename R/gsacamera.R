# Wrapper function to run CAMERA on different contrasts

runcamera <- function(voom.results, contrast, gs.annot,
                    num.workers=4, verbose = TRUE){     
    # run CAMERA and write out ranked 'gene sets' for each 'contrast'  
    if (is.matrix(contrast)){
        contr.names = colnames(contrast)
        contr.num = ncol(contrast)
    }else{
        contr.names = names(contrast)
        contr.num = length(contrast)
    }   
#    camera.results = vector("list", ncol(contrast)) 
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
#       print(args.all[[i]]) 
    }
    names(args.all) = contr.names
    if (Sys.info()['sysname'] == "Windows" || contr.num <= 1)
        camera.results = lapply(args.all, runcamera.contrast)
    else
        camera.results = mclapply(args.all, runcamera.contrast, 
                mc.cores=num.workers)
#    names(camera.results) = colnames(contrast)
    return(camera.results)
}

runcamera.contrast <- function(args){
    if (args$verbose)
        print(paste0("   Running CAMERA for ", args$contrast.name))
    else
        cat(".")
    camera.results = camera(y=args$voom.results, 
            index=args$gs.annot@idx, 
            design=args$voom.results$design, 
            contrast=args$contrast) # , allow.neg.cor=TRUE, inter.gene.cor=NA 
#       print(head(camera.results[[i]]))
    
    camera.results = 
            camera.results[order(camera.results[,"PValue"]),]
    
    camera.results = cbind(Rank=seq(1, 
                    nrow(camera.results)), camera.results)
    colnames(camera.results)[which(colnames(camera.results) == "PValue")] = "p.value"
    colnames(camera.results)[which(colnames(camera.results) == "FDR")] = "p.adj"
    return(camera.results)
}