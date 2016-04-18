# Wrapper function to run PADOG on different contrasts

runpadog <- function(voom.results, contrast, gs.annot,  ranked.gs.dir="", 
output = TRUE,
        num.workers=4){     

    # run padog and write out ranked 'gene sets' for each 'contrast'    
    file.name = paste0(ranked.gs.dir, "/padog-ranked-", gs.annot$label, 
"-gene-sets-", 
            sub(" - ", "-", colnames(contrast)), '.txt')        
    padog.results = vector("list", ncol(contrast))  
    data.log = voom.results$E
    design = voom.results$design
    rownames(data.log) = as.character(seq(1, nrow(data.log)))   
    sam.idx = 1:ncol(data.log)
    
    gsets = list()
    for (j in 1:length(gs.annot$idx)){
        gsets[[j]] = as.character(gs.annot$idx[[j]])
    }
    names(gsets) = names(gs.annot$idx)
    set.seed(05081986)
    
    for(i in 1:ncol(contrast)){
        print(paste0("   Running PADOG for ", colnames(contrast)[i]))
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
        group = c(rep("d", length(tre.sam.indx)), rep("c", 
length(cnt.sam.indx)))
        
        padog.results[[i]] = padog(esetm=data.log.sel, group=group, 
paired=FALSE, 
                gslist=gsets, NI=100, verbose=FALSE)
        
        padog.results[[i]] = 
padog.results[[i]][order(padog.results[[i]][,"Ppadog"]),]   
        padog.results[[i]] = cbind(Rank=seq(1, 
nrow(padog.results[[i]])), padog.results[[i]])
        if (!output)
            next
        writeResultsToHTML(colnames(contrast)[i], padog.results[[i]], 
gs.annot, "PADOG", file.name[i])
    }
    names(padog.results) = colnames(contrast)
    return(padog.results)
}