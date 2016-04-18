# Wrapper function to run GSVA, PLAGE, ZSCORE, and SSGSEA on different contrasts


rungsva <- function(method, voom.results, contrast, gs.annot,  
ranked.gs.dir="", output = TRUE,
        num.workers=4){     

    # run gsva and write out ranked 'gene sets' for each 'contrast'
    
    file.name = paste0(ranked.gs.dir, "/", method,"-ranked-", 
gs.annot$label, "-gene-sets-", 
            sub(" - ", "-", colnames(contrast)), '.txt')        
    gsva.results = vector("list", ncol(contrast))   
    data.log = voom.results$E
    design = voom.results$design
    rownames(data.log) = as.character(seq(1, nrow(data.log)))       
    
    gsets = list()
    for (j in 1:length(gs.annot$idx)){
        gsets[[j]] = as.character(gs.annot$idx[[j]])
    }
    names(gsets) = names(gs.annot$idx)
    
    args.all = list()
    for(i in 1:ncol(contrast)){
        args.all[[colnames(contrast)[i]]] = list(contrast=contrast,
                                        
        i = i,
                                        
        design = design,
                                        
        method=method,
                                        
        data.log = data.log,
                                        
        gsets = gsets,
                                        
        gs.annot = gs.annot,
                                        
        file.name = file.name[i],
                                        
        output = output                                                
                                        
        )
    }
    if (Sys.info()['sysname'] == "Windows" || ncol(contrast) <= 1)
        gsva.results = lapply(args.all, rungsva.contrast)
    else
        gsva.results = mclapply(args.all, rungsva.contrast, 
mc.cores=num.workers)
    #names(gsva.results) = colnames(contrast)
    return(gsva.results)
}

rungsva.contrast <- function(args){
    set.seed(519863)
    print(paste0("   Running ", toupper(args$method)," for ", 
colnames(args$contrast)[args$i]))   
    d = args$design[, args$contrast[,args$i] > 0]
    sam.idx = 1:ncol(args$data.log)
    if (is.null(ncol(d))){      
        tre.sam.indx = sam.idx[ d == 1]     
    }else if (ncol(d) > 1){     
        tre.sam.indx = c()
        for (j in 1:ncol(d))
            tre.sam.indx = c(tre.sam.indx, sam.idx[ d[,j] == 1])
    }
    else
        stop("Invalid contrasts selected.") 
    if (length(tre.sam.indx) == 1)
        tre.sam.indx = rep(tre.sam.indx, 3)
    else if (length(tre.sam.indx) < 3)
        tre.sam.indx = c(tre.sam.indx, sample(tre.sam.indx, 3 - 
length(tre.sam.indx)))
    d = args$design[, args$contrast[,args$i] < 0]
    if (is.null(ncol(d))){
        cnt.sam.indx = sam.idx[ d == 1]
    }else if (ncol(d) > 1){
        cnt.sam.indx = c()
        for (j in 1:ncol(d))
            cnt.sam.indx = c(cnt.sam.indx, sam.idx[ d[,j] == 1])
    }
    else
        stop("Invalid contrasts selected.")     
    if (length(cnt.sam.indx) == 1)
        cnt.sam.indx = rep(cnt.sam.indx, 3)
    else if (length(cnt.sam.indx) < 3)
        cnt.sam.indx = c(cnt.sam.indx, sample(cnt.sam.indx, 3 - 
length(cnt.sam.indx)))

    data.log.sel = args$data.log[, c(tre.sam.indx, cnt.sam.indx)]   

    gs.es = gsva(expr=data.log.sel, gset.idx.list=args$gsets, mx.diff=TRUE, 
min.sz=3, 
            method=args$method, parallel.sz=1, 
            verbose=FALSE, rnaseq=FALSE)#$es.obs
  
    if (args$method == "gsva"){
        gs.es = gs.es$es.obs
    }
    design.sel = args$design[c(tre.sam.indx, cnt.sam.indx), 
            args$contrast[,args$i] > 0 | args$contrast[,args$i] < 0 
]
    contrast.sel = args$contrast[args$contrast[,args$i] != 0,args$i]
    gs.fit = lmFit(gs.es, design=design.sel)
    gs.fit = contrasts.fit(gs.fit, contrast.sel)
    gs.fit = eBayes(gs.fit)
    gsva.results =  topTable(gs.fit, coef=1, , number=Inf, sort.by="p", 
adjust.method="BH") 
    gsva.results = cbind(Rank=seq(1, nrow(gsva.results)), gsva.results)
    
    
    if (args$output)
        writeResultsToHTML(colnames(args$contrast)[args$i], 
gsva.results, args$gs.annot, toupper(args$method), args$file.name)

    return(gsva.results)
}