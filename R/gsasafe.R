# Wrapper function to run SAFE (Significance Analysis of Function and 
# Expression) on different contrasts

runsafe <- function(voom.results, contrast, gs.annot,  ranked.gs.dir="", output 
= TRUE,
        num.workers=4, verbose = TRUE){     

    # run safe and write out ranked 'gene sets' for each 'contrast' 
    file.name = paste0(ranked.gs.dir, "/safe-ranked-", gs.annot$label, 
"-gene-sets-", 
            sub(" - ", "-", colnames(contrast)), '.txt')        
    safe.results = vector("list", ncol(contrast))   
    data.log = voom.results$E
    design = voom.results$design
    rownames(data.log) = as.character(seq(1, nrow(data.log)))   
    sam.idx = 1:ncol(data.log)
    
    gsets = list()
    for (j in 1:length(gs.annot$idx)){
        gsets[[j]] = as.character(gs.annot$idx[[j]])
    }
    names(gsets) = names(gs.annot$idx)  
    capture.output(C.mat <- getCmatrix(keyword.list=gsets, present.genes=rownames(data.log)))
    set.seed(05081986)
#   print(C.mat)
#   print(names(gsets)[which(!names(gsets) %in% C.mat$col.names)])
    for(i in 1:ncol(contrast)){
        if (verbose)
            print(paste0("   Running SAFE for ", colnames(contrast)[i]))
        else
            cat(".")
        d = design[, contrast[,i] > 0]
        if (is.null(ncol(d))){
            tre.sam.indx = sam.idx[ d == 1]
        }else if (ncol(d) > 1){
            tre.sam.indx = c()
            for (j in 1:ncol(d))
                tre.sam.indx = c(tre.sam.indx, sam.idx[ d[,j] 
== 1])
        }
        else
            stop("Invalid contrasts selected.")
        if (length(tre.sam.indx) == 1)
            tre.sam.indx = rep(tre.sam.indx, 3)
        else if (length(tre.sam.indx) < 3)
            tre.sam.indx = c(tre.sam.indx, sample(tre.sam.indx, 3 - 
length(tre.sam.indx)))
        d = design[, contrast[,i] < 0]
        if (is.null(ncol(d))){
            cnt.sam.indx = sam.idx[ d == 1]
        }else if (ncol(d) > 1){
            cnt.sam.indx = c()
            for (j in 1:ncol(d))
                cnt.sam.indx = c(cnt.sam.indx, sam.idx[ d[,j] 
== 1])
        }
        else
            stop("Invalid contrasts selected.")     
        if (length(cnt.sam.indx) == 1)
            cnt.sam.indx = rep(cnt.sam.indx, 3)
        else if (length(cnt.sam.indx) < 3)
            cnt.sam.indx = c(cnt.sam.indx, sample(cnt.sam.indx, 3 - 
length(cnt.sam.indx)))
        data.log.sel = data.log[, c(tre.sam.indx, cnt.sam.indx)]
        group = c(rep("1", length(tre.sam.indx)), rep("0", 
length(cnt.sam.indx)))
        
        capture.output(temp <- safe(X.mat=data.log.sel, y.vec=group,C.mat=C.mat, 
                print.it=FALSE, Pi.mat=100))#, 
# method="bootstrap.t", , method="express"
#       print(temp)
        safe.results[[i]] = safe.toptable(temp, 
number=length(gsets),description = FALSE)
        #print(head(safe.results[[i]]))
        safe.results[[i]][, "P.value"] = as.numeric(safe.results[[i]][, 
"P.value"])
        safe.results[[i]][, "Adj.p.value"] = 
as.numeric(safe.results[[i]][, "Adj.p.value"])
#       print(summary(as.numeric(safe.results[[i]][, "P.value"])))
#       stop("here")
        rownames(safe.results[[i]]) = safe.results[[i]][, "GenesetID"]
        safe.results[[i]] = safe.results[[i]][, -c(1,6)]
        safe.results[[i]] = cbind(Rank=seq(1, nrow(safe.results[[i]])), 
safe.results[[i]])
        #safe.results[[i]] = 
safe.results[[i]][order(safe.results[[i]][,"P.value"]),]    
        if (!output)
            next
        writeResultsToHTML(colnames(contrast)[i], safe.results[[i]], 
gs.annot, "safe", file.name[i])
    }
    names(safe.results) = colnames(contrast)
    return(safe.results)
}
