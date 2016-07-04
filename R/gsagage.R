# Wrapper function to run GAGE on different contrasts

rungage <- function(voom.results, contrast, gs.annot, ranked.gs.dir="", output 
= TRUE,
        num.workers=4, verbose = TRUE){     
    # run GAGE and write out ranked 'gene sets' for each 'contrast'

    file.name = paste0(ranked.gs.dir, "/gage-ranked-", gs.annot@label, 
"-gene-sets-", 
            sub(" - ", "-", colnames(contrast)), '.txt')        
    
    gage.results = vector("list", ncol(contrast))   
    data.log = voom.results$E
    design = voom.results$design
    rownames(data.log) = as.character(seq(1, nrow(data.log)))   
    sam.idx = 1:ncol(data.log)
    
    gsets = list()
    for (j in 1:length(gs.annot@idx)){
        gsets[[j]] = as.character(gs.annot@idx[[j]])
    }
    names(gsets) = names(gs.annot@idx)
    set.seed(05081986)
    #print(contrast)
    #print(design)
    for(i in 1:ncol(contrast)){
        if (verbose)
            print(paste0("   Running GAGE for ", colnames(contrast)[i]))
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
        #print(cnt.sam.indx)
        if (length(cnt.sam.indx) == 1)
            cnt.sam.indx = rep(cnt.sam.indx, 3)
        else if (length(cnt.sam.indx) < 3)
            cnt.sam.indx = c(cnt.sam.indx, sample(cnt.sam.indx, 3 - 
length(cnt.sam.indx)))
        #print(cnt.sam.indx)
        data.log.sel = data.log[, c(cnt.sam.indx, tre.sam.indx)]
        #print(head(data.log.sel))
        # same.dir=FALSE ==> Two directional test
        gage.results[[i]] = gage(exprs=data.log.sel, gsets=gsets, 
ref=seq(1,length(cnt.sam.indx)),
                samp=seq(length(cnt.sam.indx) + 
1,ncol(data.log.sel)), same.dir=FALSE, 
                compare = "unpaired")$greater[, 1:5]        
        # returns PropDown/PropUp ==> proportion of genes that are 
# down/up-regulated
        gage.results[[i]] = 
gage.results[[i]][order(gage.results[[i]][,"p.val"]),]  
        gage.results[[i]] = cbind(Rank=seq(1, nrow(gage.results[[i]])), 
gage.results[[i]])
        #print(head(gage.results[[i]]))
        
        if (!output)
            next
        writeResultsToHTML(colnames(contrast)[i], gage.results[[i]], 
gs.annot, "GAGE", file.name[i])
    }
    #stop("here")
    names(gage.results) = colnames(contrast)
    return(gage.results)
}
