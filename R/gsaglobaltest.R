# Wrapper function to run GLOBALTEST on different contrasts


runglobaltest <- function(voom.results, contrast, gs.annot,  ranked.gs.dir="", 
output = TRUE,
        num.workers=4){     

    # run globaltest and write out ranked 'gene sets' for each 'contrast'   
    file.name = paste0(ranked.gs.dir, "/globaltest-ranked-", 
gs.annot$label, "-gene-sets-", 
            sub(" - ", "-", colnames(contrast)), '.txt')        
    globaltest.results = vector("list", ncol(contrast)) 
    data.log = voom.results$E
    design = voom.results$design
    rownames(data.log) = as.character(seq(1, nrow(data.log)))   
    sam.idx = 1:ncol(data.log)
    
    gt.options(transpose=TRUE)
    
    for(i in 1:ncol(contrast)){
        print(paste0("   Running GLOBALTEST for ", 
colnames(contrast)[i]))
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
#       if (length(tre.sam.indx) < 3)
#           tre.sam.indx = c(tre.sam.indx, sample(tre.sam.indx, 3 - 
# length(tre.sam.indx)))
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
#       if (length(cnt.sam.indx) < 3)
#           cnt.sam.indx = c(cnt.sam.indx, sample(cnt.sam.indx, 3 - 
# length(cnt.sam.indx)))
        data.log.sel = data.log[, c(tre.sam.indx, cnt.sam.indx)]
        group = c(rep("d", length(tre.sam.indx)), rep("c", 
length(cnt.sam.indx)))
        colnames(data.log.sel) = group
        
        globaltest.results[[i]] = result(gt(factor(group), 
data.log.sel, subsets=gs.annot$idx, permutations=10000))
        
        globaltest.results[[i]] = 
globaltest.results[[i]][order(globaltest.results[[i]][,"p-value"]),]
        globaltest.results[[i]] = cbind(Rank=seq(1, 
nrow(globaltest.results[[i]])), globaltest.results[[i]])
        if (!output)
            next
        writeResultsToHTML(colnames(contrast)[i], 
globaltest.results[[i]], gs.annot, "GLOBALTEST", file.name[i])
    }
    names(globaltest.results) = colnames(contrast)
    return(globaltest.results)
}
