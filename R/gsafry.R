# Wrapper function to run FRY on different contrasts

runfry <- function(voom.results, contrast, gs.annot,  ranked.gs.dir="", 
        output = TRUE,
        num.workers=4, verbose = TRUE){     
    # run fry and write out ranked 'gene sets' for each 'contrast'    
    
    file.name = paste0(ranked.gs.dir, "/fry-ranked-", gs.annot$label, 
            "-gene-sets-", 
            sub(" - ", "-", colnames(contrast)), '.txt')        
    fry.results = vector("list", ncol(contrast))  
    for(i in 1:ncol(contrast)){
        if (verbose)
            print(paste0("   Running FRY for ", colnames(contrast)[i]))
        else
            cat(".")
        capture.output(fry.results[[i]] <- fry(y=voom.results, 
                index=gs.annot$idx, 
                design=voom.results$design, 
                contrast=contrast[,i]))
        # returns PropDown/PropUp ==> proportion of genes that are 
# down/up-regulated
        fry.results[[i]] = fry.results[[i]][order(fry.results[[i]][, 
                                "PValue"]),]
        fry.results[[i]] = cbind(Rank=seq(1, 
                        nrow(fry.results[[i]])), fry.results[[i]])
        if (!output)
            next
        writeResultsToHTML(colnames(contrast)[i], fry.results[[i]], 
                gs.annot, "FRY", file.name[i])
    }
    names(fry.results) = colnames(contrast)
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