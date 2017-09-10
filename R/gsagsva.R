# Wrapper function to run GSVA, PLAGE, ZSCORE, and SSGSEA on different contrasts


rungsva <- function(method, voom.results, contrast, gs.annot, 
        num.workers=4, verbose = TRUE){     

    # run gsva and write out ranked 'gene sets' for each 'contrast'
    if (is.matrix(contrast)){
        contr.names = colnames(contrast)
        contr.num = ncol(contrast)
    }else{
        contr.names = names(contrast)
        contr.num = length(contrast)
    }
    
    gsets = list()
    for (j in 1:length(gs.annot@idx)){
        gsets[[j]] = as.character(gs.annot@idx[[j]])
    }
    names(gsets) = names(gs.annot@idx)

    if (verbose)
        message("   Calculating gene set-level stats using ", 
                        toupper(method))
    else
        message(".", appendLF = FALSE)
    # transform scores in gene set space using parallel computing
    data.log = voom.results$E
    rownames(data.log) = as.character(seq(1, nrow(data.log)))
    gs.es = calculateSetScores.parallel(data.log, gsets, method, num.workers)
    # fit the gene set scores and find DE gene sets
    gs.fit = lmFit(gs.es, design=voom.results$design)
    if (is.matrix(contrast)){
        gs.fit = contrasts.fit(gs.fit, contrast)
        coefs = 1:contr.num
    }else{
        coefs = contrast
    }
    gs.fit = eBayes(gs.fit)
    gsva.results = vector("list", contr.num)
    for(i in 1:contr.num){
        if (verbose)
            message("   Running ", toupper(method)," for ", 
                            contr.names[i])
        else
            message(".", appendLF = FALSE)        
        gsva.results[[i]] =  topTable(gs.fit, coef=coefs[i], , number=Inf, sort.by="p", 
                adjust.method="BH") 
        gsva.results[[i]] = gsva.results[[i]][order(gsva.results[[i]][, "P.Value"], 
                        -gsva.results[[i]][, "B"]),]
        gsva.results[[i]] = cbind(Rank=seq(1, nrow(gsva.results[[i]])), 
                gsva.results[[i]])        
        
        colnames(gsva.results[[i]])[which(
                        colnames(gsva.results[[i]]) == "P.Value")] = "p.value"
        colnames(gsva.results[[i]])[which(
                        colnames(gsva.results[[i]]) == "adj.P.Val")] = "p.adj"
    }
    names(gsva.results) = contr.names
    return(gsva.results)
}

calculateSetScores.parallel <- function(data.log, gsets, method, num.workers){
    args.all = list()
    sets.per.task = 50
    total.tasks = ceiling(length(gsets) /  sets.per.task)    
    for (i in 1:total.tasks){
        gsetsi = gsets[((i-1)*sets.per.task + 1):(i*sets.per.task)]
        gsetsi = gsetsi[!sapply(gsetsi, is.null)]
        args.all[[paste0("task", i)]] = list(
                data.log = data.log,
                gsets = gsetsi,
                method = method
        )
    }
    # parallelize the calculation of gene set scores to speed up the algorithms
    if (Sys.info()['sysname'] == "Windows")
        gs.es.all = lapply(args.all, rungsva.subcollection)
    else
        gs.es.all = mclapply(args.all, rungsva.subcollection, 
                mc.cores=num.workers)
    # collect gene set scores from differet workers
    gs.es = c()
    for (i in 1:total.tasks)
        gs.es = rbind(gs.es, gs.es.all[[paste0("task", i)]])
    rownames(gs.es) = names(gsets) 
    return(gs.es)
}

rungsva.subcollection <- function(args){
    set.seed(519863)
    gs.es = gsva(expr=args$data.log, gset.idx.list=args$gsets, mx.diff=TRUE, 
            min.sz=1, 
            method=args$method, parallel.sz=1, 
            verbose=FALSE, kcdf = "Gaussian")#$es.obs

    # Depreciated in Version 1.25.6  
    # if (args$method == "gsva"){
    #     gs.es = gs.es$es.obs
    # }
    
    return(gs.es)
}








